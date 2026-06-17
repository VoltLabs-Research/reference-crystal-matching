#include "site_assignment.h"
#include "statistics.h"

#include <volt/structures/crystal_structure_types.h>

#include <spdlog/spdlog.h>

#include <algorithm>
#include <cmath>
#include <unordered_map>

namespace Volt{

namespace{

constexpr int kNeighborSlots = MAX_NEIGHBORS;

// The central 20%-80% box of the frame, in lab coordinates. Both the snap-
// residual and loop-residual phases restrict sampling to this interior to keep
// surface/boundary atoms out of the quality metrics.
struct BulkRegion{
    Point3 boundsLow;
    double spanX = 0, spanY = 0, spanZ = 0;

    static BulkRegion of(const LammpsParser::Frame& frame, std::size_t atomCount){
        BulkRegion region;
        const AffineTransformation& cellMatrix = frame.simulationCell.matrix();
        region.spanX = cellMatrix.column(0).length();
        region.spanY = cellMatrix.column(1).length();
        region.spanZ = cellMatrix.column(2).length();
        region.boundsLow = frame.positions[0];
        for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
            for(int axis = 0; axis < 3; ++axis){
                region.boundsLow[axis] = std::min(region.boundsLow[axis], frame.positions[atomIndex][axis]);
            }
        }
        return region;
    }

    bool contains(const Point3& position) const{
        return !(position.x() < boundsLow.x() + 0.2 * spanX || position.x() > boundsLow.x() + 0.8 * spanX ||
                 position.y() < boundsLow.y() + 0.2 * spanY || position.y() > boundsLow.y() + 0.8 * spanY ||
                 position.z() < boundsLow.z() + 0.2 * spanZ || position.z() > boundsLow.z() + 0.8 * spanZ);
    }
};

// Local slot of `targetAtom` within `fromAtom`'s neighbour list, or -1.
int slotInto(const AllAtomNeighbors& neighbors, std::size_t fromAtom, int targetAtom){
    const int start = neighbors.offsets[fromAtom];
    const int end = neighbors.offsets[fromAtom + 1];
    for(int slot = start; slot < end; ++slot){
        if(neighbors.indices[static_cast<std::size_t>(slot)] == targetAtom){
            return slot - start;
        }
    }
    return -1;
}

Vector3& overrideAt(std::vector<Vector3>& overrides, std::size_t atomIndex, int slot){
    return overrides[atomIndex * static_cast<std::size_t>(kNeighborSlots) + static_cast<std::size_t>(slot)];
}

// Phase 1 — assign every atom to the reference basis site whose neighbour shell
// (in the grain frame) best matches the atom's measured neighbourhood.
std::vector<int> assignBasisSites(
    const SiteMatchInputs& inputs,
    double cutoff,
    const AllAtomNeighbors& neighbors,
    AssignmentStats& stats)
{
    const std::size_t atomCount = static_cast<std::size_t>(inputs.frame.natoms);
    const PerfectReference& reference = inputs.reference;

    std::unordered_map<int, std::vector<int>> sitesBySpecies;
    for(int siteIndex = 0; siteIndex < static_cast<int>(reference.sites.size()); ++siteIndex){
        sitesBySpecies[reference.sites[static_cast<std::size_t>(siteIndex)].species].push_back(siteIndex);
    }

    stats.minNeighborCount = std::numeric_limits<int>::max();
    stats.maxNeighborCount = 0;
    stats.atomsWithZeroNeighbors = 0;
    stats.edgesTotal = 0;
    stats.atomsBySpeciesUnassigned = 0;

    std::vector<int> basisSiteOfAtom(atomCount, -1);
    for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
        const int start = neighbors.offsets[atomIndex];
        const int end = neighbors.offsets[atomIndex + 1];
        const int neighborCount = end - start;
        stats.minNeighborCount = std::min(stats.minNeighborCount, neighborCount);
        stats.maxNeighborCount = std::max(stats.maxNeighborCount, neighborCount);
        if(neighborCount == 0){
            ++stats.atomsWithZeroNeighbors;
            continue;
        }

        const int centralSpecies = inputs.atomSpecies[atomIndex];
        std::vector<Vector3> crystalFrameDeltas(static_cast<std::size_t>(neighborCount));
        for(int slot = 0; slot < neighborCount; ++slot){
            crystalFrameDeltas[static_cast<std::size_t>(slot)] =
                inputs.grainRotationTransposed * neighbors.deltas[static_cast<std::size_t>(start + slot)];
        }

        int bestSite = -1;
        double bestScore = std::numeric_limits<double>::max();
        const auto speciesSites = sitesBySpecies.find(centralSpecies);
        if(speciesSites != sitesBySpecies.end()){
            for(int siteIndex : speciesSites->second){
                const BasisSite& site = reference.sites[static_cast<std::size_t>(siteIndex)];
                double score = 0;
                for(int slot = 0; slot < neighborCount; ++slot){
                    const int neighborSpecies =
                        inputs.atomSpecies[static_cast<std::size_t>(neighbors.indices[static_cast<std::size_t>(start + slot)])];
                    double bestShellDistance = std::numeric_limits<double>::max();
                    for(const auto& [shellSpecies, shellVector] : site.shell){
                        if(shellSpecies != neighborSpecies){
                            continue;
                        }
                        bestShellDistance = std::min(bestShellDistance,
                            (shellVector - crystalFrameDeltas[static_cast<std::size_t>(slot)]).squaredLength());
                    }
                    if(bestShellDistance < std::numeric_limits<double>::max()){
                        score += bestShellDistance;
                    }else{
                        score += cutoff * cutoff;
                    }
                }
                if(score < bestScore){
                    bestScore = score;
                    bestSite = siteIndex;
                }
            }
        }

        stats.edgesTotal += neighborCount;
        if(bestSite < 0){
            ++stats.atomsBySpeciesUnassigned;
            continue;
        }
        basisSiteOfAtom[atomIndex] = bestSite;
    }
    if(stats.minNeighborCount == std::numeric_limits<int>::max()){
        stats.minNeighborCount = 0;
    }
    return basisSiteOfAtom;
}

// Phase 2 — for every edge between two assigned atoms, derive the ideal lab-frame
// neighbour vector from the reference basis (nearest periodic image) and store it
// as the override. Returns bulk snap residuals (measured vs ideal length).
std::vector<double> computeIdealOverrides(
    const SiteMatchInputs& inputs,
    const AllAtomNeighbors& neighbors,
    const std::vector<int>& basisSiteOfAtom,
    const BulkRegion& bulk,
    std::vector<Vector3>& overrides,
    AssignmentStats& stats)
{
    const std::size_t atomCount = static_cast<std::size_t>(inputs.frame.natoms);
    const PerfectReference& reference = inputs.reference;
    const double cellLengths[3] = {reference.cellLengthA, reference.cellLengthB, reference.cellLengthC};

    stats.edgesWithIdealVector = 0;
    std::vector<double> snapResiduals;

    for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
        const int basisSite = basisSiteOfAtom[atomIndex];
        if(basisSite < 0){
            continue;
        }
        const Vector3& atomFractional = reference.sites[static_cast<std::size_t>(basisSite)].fractional;
        const int start = neighbors.offsets[atomIndex];
        const int end = neighbors.offsets[atomIndex + 1];
        const bool atomInBulk = bulk.contains(inputs.frame.positions[atomIndex]);

        for(int slot = start; slot < end; ++slot){
            const int neighbor = neighbors.indices[static_cast<std::size_t>(slot)];
            const int neighborBasisSite = basisSiteOfAtom[static_cast<std::size_t>(neighbor)];
            if(neighborBasisSite < 0){
                continue;
            }
            const Vector3& neighborFractional = reference.sites[static_cast<std::size_t>(neighborBasisSite)].fractional;
            const Vector3 crystalFrameDelta = inputs.grainRotationTransposed * neighbors.deltas[static_cast<std::size_t>(slot)];

            double idealFractionalDelta[3];
            for(int axis = 0; axis < 3; ++axis){
                const double measuredFractional = crystalFrameDelta[axis] / cellLengths[axis];
                const double basisFractional = neighborFractional[axis] - atomFractional[axis];
                const double imageOffset = std::round(measuredFractional - basisFractional);
                idealFractionalDelta[axis] = (basisFractional + imageOffset) * cellLengths[axis];
            }
            const Vector3 idealLabFrame = inputs.grainRotation *
                Vector3(idealFractionalDelta[0], idealFractionalDelta[1], idealFractionalDelta[2]);

            overrideAt(overrides, atomIndex, slot - start) = idealLabFrame;
            ++stats.edgesWithIdealVector;
            if(atomInBulk && snapResiduals.size() < 20000){
                snapResiduals.push_back((neighbors.deltas[static_cast<std::size_t>(slot)] - idealLabFrame).length());
            }
        }
    }
    return snapResiduals;
}

// Phase 3 — force the elastic mapping to be reciprocal (v_AB = -v_BA) on every
// bidirectional edge, copying whichever direction was assigned onto its reverse.
// This is the invariant Burgers circuits depend on, so it is self-checked.
void enforceReciprocity(
    const AllAtomNeighbors& neighbors,
    std::vector<Vector3>& overrides,
    AssignmentStats& stats)
{
    const std::size_t atomCount = neighbors.offsets.size() - 1;

    for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
        const int start = neighbors.offsets[atomIndex];
        const int end = neighbors.offsets[atomIndex + 1];
        for(int slot = start; slot < end; ++slot){
            const int neighbor = neighbors.indices[static_cast<std::size_t>(slot)];
            if(neighbor <= static_cast<int>(atomIndex)){
                continue;
            }
            const int slotForward = slot - start;
            const int slotReverse = slotInto(neighbors, static_cast<std::size_t>(neighbor), static_cast<int>(atomIndex));
            if(slotReverse < 0){
                continue;
            }
            Vector3& forwardVector = overrideAt(overrides, atomIndex, slotForward);
            Vector3& reverseVector = overrideAt(overrides, static_cast<std::size_t>(neighbor), slotReverse);
            const bool forwardZero = forwardVector.isZero(1e-9);
            const bool reverseZero = reverseVector.isZero(1e-9);
            if(!forwardZero){
                reverseVector = -forwardVector;
            }else if(!reverseZero){
                forwardVector = -reverseVector;
            }
        }
    }

    // Self-check: recount assigned edges and verify reciprocity on every
    // bidirectional edge. asserts compile out under NDEBUG (this ships Release),
    // so a surviving violation is logged loudly instead.
    stats.edgesWithIdealVector = 0;
    long reciprocityViolations = 0;
    for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
        const int start = neighbors.offsets[atomIndex];
        const int end = neighbors.offsets[atomIndex + 1];
        for(int slot = start; slot < end; ++slot){
            const Vector3& forwardVector = overrideAt(overrides, atomIndex, slot - start);
            if(forwardVector.isZero(1e-9)){
                continue;
            }
            ++stats.edgesWithIdealVector;
            const int neighbor = neighbors.indices[static_cast<std::size_t>(slot)];
            const int slotReverse = slotInto(neighbors, static_cast<std::size_t>(neighbor), static_cast<int>(atomIndex));
            if(slotReverse < 0){
                continue;
            }
            const Vector3& reverseVector = overrideAt(overrides, static_cast<std::size_t>(neighbor), slotReverse);
            if(!reverseVector.isZero(1e-9) && !(forwardVector + reverseVector).isZero(1e-6)){
                ++reciprocityViolations;
            }
        }
    }
    if(reciprocityViolations > 0){
        spdlog::error("FullCrystalContext: reciprocity self-check FAILED on {} edges (v_AB != -v_BA); "
                      "Burgers circuits will not close", reciprocityViolations);
    }
}

// Phase 4 — sample closed 3-edge loops (atom -> n1 -> n2 -> atom) in the bulk and
// record their ideal-vector closure error; a measure of how well the mapping
// satisfies the discrete compatibility condition.
void sampleLoopResiduals(
    const SiteMatchInputs& inputs,
    const AllAtomNeighbors& neighbors,
    const std::vector<Vector3>& overrides,
    const BulkRegion& bulk,
    AssignmentStats& stats)
{
    const std::size_t atomCount = static_cast<std::size_t>(inputs.frame.natoms);
    const auto overrideOf = [&](std::size_t atomIndex, int slot) -> Vector3{
        return overrides[atomIndex * static_cast<std::size_t>(kNeighborSlots) + static_cast<std::size_t>(slot)];
    };

    std::vector<double> loopResiduals;
    for(std::size_t atomIndex = 0; atomIndex < atomCount && loopResiduals.size() < 2000; ++atomIndex){
        if(!bulk.contains(inputs.frame.positions[atomIndex])){
            continue;
        }
        const int start = neighbors.offsets[atomIndex];
        const int end = neighbors.offsets[atomIndex + 1];
        for(int firstSlot = start; firstSlot < end; ++firstSlot){
            const int firstNeighbor = neighbors.indices[static_cast<std::size_t>(firstSlot)];
            for(int secondSlot = firstSlot + 1; secondSlot < end; ++secondSlot){
                const int secondNeighbor = neighbors.indices[static_cast<std::size_t>(secondSlot)];
                const int slotFirstToSecond = slotInto(neighbors, static_cast<std::size_t>(firstNeighbor), secondNeighbor);
                if(slotFirstToSecond < 0){
                    continue;
                }
                const int slotSecondToAtom = slotInto(neighbors, static_cast<std::size_t>(secondNeighbor), static_cast<int>(atomIndex));
                if(slotSecondToAtom < 0){
                    continue;
                }
                const Vector3 edgeAtomToFirst = overrideOf(atomIndex, firstSlot - start);
                const Vector3 edgeFirstToSecond = overrideOf(static_cast<std::size_t>(firstNeighbor), slotFirstToSecond);
                const Vector3 edgeSecondToAtom = overrideOf(static_cast<std::size_t>(secondNeighbor), slotSecondToAtom);
                if(edgeAtomToFirst.isZero(1e-9) || edgeFirstToSecond.isZero(1e-9) || edgeSecondToAtom.isZero(1e-9)){
                    continue;
                }
                loopResiduals.push_back((edgeAtomToFirst + edgeFirstToSecond + edgeSecondToAtom).length());
                if(loopResiduals.size() >= 2000){
                    break;
                }
            }
            if(loopResiduals.size() >= 2000){
                break;
            }
        }
    }
    stats.bulkLoopsSampled = static_cast<int>(loopResiduals.size());
    stats.bulkLoopMedianResidual = loopResiduals.empty() ? -1 : medianOf(loopResiduals);
}

} // namespace

CutoffAssignment assignForCutoff(const SiteMatchInputs& inputs, double cutoff){
    const std::size_t atomCount = static_cast<std::size_t>(inputs.frame.natoms);

    CutoffAssignment out;
    out.neighbors = buildAllAtomNeighbors(inputs.frame.positions.data(), atomCount, inputs.frame.simulationCell, cutoff);
    out.overrides.assign(atomCount * static_cast<std::size_t>(kNeighborSlots), Vector3::Zero());

    const BulkRegion bulk = BulkRegion::of(inputs.frame, atomCount);

    out.basisSiteOfAtom = assignBasisSites(inputs, cutoff, out.neighbors, out.stats);
    const std::vector<double> snapResiduals =
        computeIdealOverrides(inputs, out.neighbors, out.basisSiteOfAtom, bulk, out.overrides, out.stats);
    enforceReciprocity(out.neighbors, out.overrides, out.stats);
    sampleLoopResiduals(inputs, out.neighbors, out.overrides, bulk, out.stats);

    out.snapResidualP99 = snapResiduals.empty()
        ? std::numeric_limits<double>::max()
        : percentileOf(snapResiduals, 0.99);
    return out;
}

}
