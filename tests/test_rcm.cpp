// Reference Crystal Matching — generic test suite.
//
// Proves RCM is NOT forsterite-specific: it registers four distinct crystals
// (cubic/orthorhombic, 2/3 species, isotropic/anisotropic) against a perfect
// reference and produces a reciprocal elastic mapping. Also locks the
// modularization — especially the v_AB = -v_BA invariant the Burgers circuits
// depend on, which the in-library self-check only logs but never enforces.
//
// IMPORTANT: this ships in a Release/-DNDEBUG build, so assert() is a no-op.
// Every check goes through CHECK(), which counts failures; main() returns
// non-zero if any fired. An assert-based test here would pass unconditionally.
//
// Fixtures under tests/fixtures/ are committed (atomsk is NOT a build/CI dep).
// To regenerate:
//   atomsk --create rocksalt   4.21  Mg O    mgo_ref.lmp
//   atomsk --create fluorite   5.46  Ca F    caf2_ref.lmp
//   atomsk --create perovskite 3.905 Sr Ti O srtio3_ref.lmp
//   cp ../../structures/cementite_ref.lmp     cementite_ref.lmp
//   atomsk mgo_ref.lmp       -duplicate 6 6 6 -rotate com z 27 -rotate com x 13 mgo_super.lmp
//   atomsk caf2_ref.lmp      -duplicate 6 6 6 -rotate com z 31                   caf2_super.lmp
//   atomsk srtio3_ref.lmp    -duplicate 6 6 6 -rotate com y 19                   srtio3_super.lmp
//   atomsk cementite_ref.lmp -duplicate 5 4 6 -rotate com z 23                   cementite_super.lmp

// The library compiles every TU with <volt/core/volt.h> force-included via a
// PRIVATE precompiled header; that PCH does not reach this test target, and the
// math headers assume volt.h was seen first. Include it before any project
// header to replicate that order (otherwise Matrix3 resolves incomplete).
#include <volt/core/volt.h>

#include "grain_frame.h"
#include "reference_lattice.h"
#include "site_assignment.h"

#include <volt/analysis/reference_crystal_matching_service.h>
#include <volt/core/lammps_parser.h>
#include <volt/structures/crystal_structure_types.h>

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

using namespace Volt;

namespace{

int g_failures = 0;

#define CHECK(cond, msg) do{ \
    if(!(cond)){ \
        std::fprintf(stderr, "FAIL %s:%d [%s]: %s\n", __FILE__, __LINE__, fx.name, msg); \
        ++g_failures; \
    } \
}while(0)

constexpr int kSlots = MAX_NEIGHBORS;

struct Fixture{
    const char* name;
    const char* refFile;
    const char* superFile;
    int anchorSpecies;
    int expectedBasisSites;   // atom count of the conventional reference cell
    bool isotropic;           // cubic -> metric_rescale (1,1,1)
};

const Fixture kFixtures[] = {
    {"rocksalt-MgO",    "mgo_ref.lmp",       "mgo_super.lmp",       1,  8, true },
    {"fluorite-CaF2",   "caf2_ref.lmp",      "caf2_super.lmp",      1, 12, true },
    {"perovskite-SrTiO3","srtio3_ref.lmp",   "srtio3_super.lmp",    1,  5, true },
    {"cementite-Fe3C",  "cementite_ref.lmp", "cementite_super.lmp", 1, 16, false},
};

bool loadFrame(const std::string& path, LammpsParser::Frame& frame){
    LammpsParser parser;
    return parser.parseFile(path, frame);
}

// Build the anchor sublattice from a parsed reference, mirroring the logic in
// full_crystal_context.cpp (filter by anchor species, dedup fractional coords).
AnchorReference anchorFrom(const PerfectReference& ref, int anchorSpecies){
    AnchorReference anchor;
    anchor.cellLengthA = ref.cellLengthA;
    anchor.cellLengthB = ref.cellLengthB;
    anchor.cellLengthC = ref.cellLengthC;
    anchor.targetSpecies = anchorSpecies;
    for(const BasisSite& site : ref.sites){
        if(site.species != anchorSpecies){
            continue;
        }
        Vector3 fractional = site.fractional;
        for(int axis = 0; axis < 3; ++axis){
            fractional[axis] -= std::floor(fractional[axis]);
        }
        bool duplicate = false;
        for(const Vector3& existing : anchor.idealFractional){
            if((existing - fractional).length() < 0.02){
                duplicate = true;
                break;
            }
        }
        if(!duplicate){
            anchor.idealFractional.push_back(fractional);
        }
    }
    return anchor;
}

void runFixture(const std::string& dir, const Fixture& fx){
    const std::string refPath   = dir + "/" + fx.refFile;
    const std::string superPath = dir + "/" + fx.superFile;

    // (a) reference_lattice: parse + cell geometry.
    const PerfectReference ref = buildPerfectReference(refPath, 0.9 * 4.0);
    CHECK(ref.ok, "buildPerfectReference failed");
    if(!ref.ok){
        return; // nothing else is meaningful without a reference
    }
    CHECK(static_cast<int>(ref.sites.size()) == fx.expectedBasisSites,
          "unexpected basis site count");
    CHECK(ref.cellLengthA > 1e-3 && ref.cellLengthB > 1e-3 && ref.cellLengthC > 1e-3,
          "degenerate reference cell lengths");
    if(fx.isotropic){
        CHECK(std::abs(ref.cellLengthA - ref.cellLengthB) < 1e-6 &&
              std::abs(ref.cellLengthA - ref.cellLengthC) < 1e-6,
              "cubic reference cell not equal-sided");
    }

    const AnchorReference anchor = anchorFrom(ref, fx.anchorSpecies);
    CHECK(!anchor.idealFractional.empty(), "no anchor sites in reference");

    // (b) grain_frame: recover orientation from the rotated supercell.
    LammpsParser::Frame superFrame;
    CHECK(loadFrame(superPath, superFrame), "cannot open supercell fixture");
    if(superFrame.natoms <= 0){
        std::fprintf(stderr, "FAIL [%s]: supercell has no atoms\n", fx.name);
        ++g_failures;
        return;
    }

    // .lmp carries species in frame.types (column "type").
    const std::vector<int>& species = superFrame.types;
    long anchorAtomsExpected = 0;
    for(int t : species){
        if(t == fx.anchorSpecies){ ++anchorAtomsExpected; }
    }
    CHECK(anchorAtomsExpected > 50, "fixture has too few anchor atoms (need >50)");

    const GrainFrame gf = recoverGrainFrame(
        superFrame.positions.data(), species.data(),
        static_cast<std::size_t>(superFrame.natoms), superFrame.simulationCell, anchor);
    CHECK(gf.ok, "recoverGrainFrame failed");
    if(!gf.ok){
        return;
    }
    CHECK(gf.bulkSnapResidual < 0.2, "grain residual too high (>0.2 A)");
    CHECK(gf.rotation.isOrthogonalMatrix(1e-6), "recovered frame not orthogonal");
    CHECK(std::abs(gf.rotation.determinant() - 1.0) < 1e-6, "recovered frame not a proper rotation");

    // (d) anchor count consistency (no magic literal).
    CHECK(gf.anchorAtomCount == static_cast<int>(anchorAtomsExpected),
          "grain frame anchor count != counted anchor atoms");

    // (c) site_assignment: reciprocity v_AB = -v_BA (THE core invariant).
    const SiteMatchInputs inputs{
        superFrame, ref, species, gf.rotation, gf.rotation.transposed()};
    const double cutoff = 0.85 * std::min({ref.cellLengthA, ref.cellLengthB, ref.cellLengthC});
    const CutoffAssignment asg = assignForCutoff(inputs, cutoff);

    CHECK(asg.snapResidualP99 < 1e30, "no cutoff produced a finite snap residual");
    CHECK(asg.stats.edgesWithIdealVector > 0, "no ideal vectors assigned");
    CHECK(asg.basisSiteOfAtom.size() == static_cast<std::size_t>(superFrame.natoms),
          "basisSiteOfAtom size mismatch");

    const auto overrideAt = [&](std::size_t atom, int slot) -> const Vector3&{
        return asg.overrides[atom * static_cast<std::size_t>(kSlots) + static_cast<std::size_t>(slot)];
    };
    const auto slotInto = [&](std::size_t fromAtom, int targetAtom) -> int{
        const int start = asg.neighbors.offsets[fromAtom];
        const int end = asg.neighbors.offsets[fromAtom + 1];
        for(int s = start; s < end; ++s){
            if(asg.neighbors.indices[static_cast<std::size_t>(s)] == targetAtom){
                return s - start;
            }
        }
        return -1;
    };

    long reciprocityViolations = 0;
    long bidirectionalEdges = 0;
    const std::size_t atomCount = static_cast<std::size_t>(superFrame.natoms);
    for(std::size_t atom = 0; atom < atomCount; ++atom){
        const int start = asg.neighbors.offsets[atom];
        const int end = asg.neighbors.offsets[atom + 1];
        for(int s = start; s < end; ++s){
            const int neighbor = asg.neighbors.indices[static_cast<std::size_t>(s)];
            if(neighbor <= static_cast<int>(atom)){
                continue;
            }
            const int slotFwd = s - start;
            const int slotRev = slotInto(static_cast<std::size_t>(neighbor), static_cast<int>(atom));
            if(slotRev < 0){
                continue;
            }
            const Vector3& fwd = overrideAt(atom, slotFwd);
            const Vector3& rev = overrideAt(static_cast<std::size_t>(neighbor), slotRev);
            if(fwd.isZero(1e-9) && rev.isZero(1e-9)){
                continue;
            }
            ++bidirectionalEdges;
            if(!(fwd + rev).isZero(1e-6)){
                ++reciprocityViolations;
            }
        }
    }
    CHECK(bidirectionalEdges > 0, "no bidirectional ideal edges to check reciprocity");
    CHECK(reciprocityViolations == 0, "reciprocity v_AB = -v_BA violated on some edges");

    // (e) full service path: metric_rescale isotropy.
    ReferenceCrystalMatchingService service;
    service.setReferenceCrystal(refPath);
    service.setAnchorSpecies(fx.anchorSpecies);
    service.setSpeciesColumn("type"); // .lmp uses the "type" column
    service.setTopologyName(fx.name);
    const std::string outBase =
        (std::filesystem::temp_directory_path() / (std::string("rcm_test_") + fx.name)).string();
    const json result = service.compute(superFrame, outBase, superPath);
    CHECK(result.contains("is_failed") && result["is_failed"] == false,
          "service.compute reported failure");
    if(result.contains("metric_rescale")){
        const double sx = result["metric_rescale"]["sx"].get<double>();
        const double sy = result["metric_rescale"]["sy"].get<double>();
        const double sz = result["metric_rescale"]["sz"].get<double>();
        if(fx.isotropic){
            CHECK(std::abs(sx - 1.0) < 1e-6 && std::abs(sy - 1.0) < 1e-6 && std::abs(sz - 1.0) < 1e-6,
                  "cubic crystal did not yield isotropic metric_rescale (1,1,1)");
        }else{
            // negative control: anisotropic cell must NOT be (1,1,1).
            CHECK(!(std::abs(sx - 1.0) < 1e-6 && std::abs(sy - 1.0) < 1e-6 && std::abs(sz - 1.0) < 1e-6),
                  "anisotropic crystal wrongly yielded isotropic metric_rescale");
        }
    }
}

} // namespace

int main(int argc, char** argv){
    const std::string dir = (argc > 1) ? argv[1] : RCM_FIXTURE_DIR;
    std::fprintf(stderr, "RCM test suite — fixtures: %s\n", dir.c_str());

    for(const Fixture& fx : kFixtures){
        runFixture(dir, fx);
    }

    if(g_failures == 0){
        std::fprintf(stderr, "OK — all RCM fixtures passed\n");
        return 0;
    }
    std::fprintf(stderr, "%d CHECK(s) FAILED\n", g_failures);
    return 1;
}
