#include <volt/analysis/reference_crystal_matching_log.h>

#include <volt/analysis/crystal_context.h>
#include <volt/registration/grain_frame.h>
#include <volt/registration/grain_segmentation.h>
#include <volt/registration/reference_lattice.h>
#include <volt/registration/site_assignment.h>

#include <spdlog/spdlog.h>

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <map>
#include <numbers>
#include <utility>

namespace Volt{

std::string ReferenceCrystalMatchingLog::describeSpeciesCounts(const std::vector<int>& atomSpecies){
    constexpr std::size_t kMaxListedSpecies = 8;
    std::vector<std::pair<int, std::size_t>> counts;
    std::size_t unlistedAtoms = 0;
    for(int species : atomSpecies){
        auto entry = std::find_if(counts.begin(), counts.end(), [species](const auto& item){
            return item.first == species;
        });
        if(entry != counts.end()){
            ++entry->second;
        }else if(counts.size() < kMaxListedSpecies){
            counts.emplace_back(species, 1);
        }else{
            ++unlistedAtoms;
        }
    }
    std::sort(counts.begin(), counts.end());

    std::string text;
    char item[64];
    for(const auto& [species, count] : counts){
        std::snprintf(item, sizeof(item), "%s%d: %zu (%.1f%%)", text.empty() ? "" : ", ",
                      species, count, percentOf(count, atomSpecies.size()));
        text += item;
    }
    if(unlistedAtoms > 0){
        std::snprintf(item, sizeof(item), ", other species: %zu (%.1f%%)",
                      unlistedAtoms, percentOf(unlistedAtoms, atomSpecies.size()));
        text += item;
    }
    return text;
}

std::string ReferenceCrystalMatchingLog::describeSiteSpecies(const std::vector<BasisSite>& sites){
    std::map<int, int> sitesPerSpecies;
    for(const BasisSite& site : sites){
        ++sitesPerSpecies[site.species];
    }
    std::string text;
    for(const auto& [species, count] : sitesPerSpecies){
        text += (text.empty() ? "species " : ", ") + std::to_string(species) + ": " + std::to_string(count);
    }
    return text;
}

std::string ReferenceCrystalMatchingLog::describeGrainSizes(const std::vector<int>& grainOfAtom, int grainCount){
    constexpr int kMaxListedGrainSizes = 8;
    std::vector<std::size_t> atomsPerGrain(static_cast<std::size_t>(std::max(grainCount, 0)), 0);
    for(int grain : grainOfAtom){
        if(grain >= 0 && grain < grainCount){
            ++atomsPerGrain[static_cast<std::size_t>(grain)];
        }
    }
    std::string text;
    char item[48];
    for(int grain = 0; grain < std::min(grainCount, kMaxListedGrainSizes); ++grain){
        const std::size_t count = atomsPerGrain[static_cast<std::size_t>(grain)];
        std::snprintf(item, sizeof(item), "%s%zu (%.1f%%)", text.empty() ? "" : ", ",
                      count, percentOf(count, grainOfAtom.size()));
        text += item;
    }
    if(grainCount > kMaxListedGrainSizes){
        text += ", ...";
    }
    return text;
}

void ReferenceCrystalMatchingLog::printReciprocityViolations(const AssignmentStats& stats){
    if(stats.reciprocityViolations > 0){
        spdlog::error("Reciprocity self-check FAILED on {} edges (v_AB != -v_BA); "
                      "Burgers circuits will not close", stats.reciprocityViolations);
    }
}

void ReferenceCrystalMatchingLog::dump(const ReferenceCrystalMatchingService::ExoticReference& reference){
    spdlog::info("Reference topology '{}' resolved to {} (anchor species {})",
                 reference.name,
                 std::filesystem::path(reference.referenceCrystal).filename().string(),
                 reference.anchorSpecies);
}

void ReferenceCrystalMatchingLog::dump(const std::vector<int>& atomSpecies){
    spdlog::info("Atom species: {}", describeSpeciesCounts(atomSpecies));
}

void ReferenceCrystalMatchingLog::dump(const PerfectReference& reference){
    spdlog::info("Reference crystal: cell {:.3f} x {:.3f} x {:.3f}, {} basis sites ({}), shell reach {:.2f} A",
                 reference.cellLengthA, reference.cellLengthB, reference.cellLengthC,
                 reference.sites.size(), describeSiteSpecies(reference.sites), reference.shellReach);
}

void ReferenceCrystalMatchingLog::dump(const GrainFrame& grainFrame){
    spdlog::info("Grain frame recovered in {} from {} anchors: rotation {:.2f} deg, median bond residual {:.3f} A",
                 phaseElapsed(), grainFrame.anchorAtomCount,
                 grainFrame.rotation.rotationAngle() * 180.0 / std::numbers::pi, grainFrame.bulkSnapResidual);
}

void ReferenceCrystalMatchingLog::dump(const GrainSegmentation& grains){
    if(grains.anchorSearch){
        const GrainSegmentation::AnchorSearch& search = *grains.anchorSearch;
        const std::size_t listedGrains = std::min(search.grainAnchors.size(), kMaxListedGrains);
        for(std::size_t grain = 0; grain < listedGrains; ++grain){
            const std::size_t anchors = search.grainAnchors[grain];
            spdlog::info("  Grain {}: {} anchors ({:.1f}%)", grain + 1, anchors, percentOf(anchors, search.anchors));
        }
        if(search.grainAnchors.empty()){
            spdlog::warn("No anchor fits the global grain frame; treating the sample as a single grain");
        }else{
            spdlog::info("  Grain search stopped with {} anchors unassigned: {}",
                         search.unassignedAnchors, search.stopReason);
        }
    }
    if(grains.unreachedAtoms > 0){
        spdlog::warn("{} atoms are not connected to any grain within {:.1f} A; assigned to grain 1",
                     grains.unreachedAtoms, AnchorReference::physicalBondCutoff);
    }
    spdlog::info("Segmented {} grain(s) in {}; atoms per grain: {}",
                 grains.grainCount, phaseElapsed(), describeGrainSizes(grains.grainOfAtom, grains.grainCount));
    _cutoffSelection.restart();
}

void ReferenceCrystalMatchingLog::dump(const CutoffAssignment& trial){
    printReciprocityViolations(trial.stats);
    spdlog::info("  Cutoff {:.2f} A: snap p90 {:.3f} A, p99 {:.3f} A, max coordination {}, {} ideal edges ({})",
                 trial.cutoff, trial.snapResidualP90, trial.snapResidualP99,
                 trial.stats.maxNeighborCount, trial.stats.edgesWithIdealVector, phaseElapsed());
}

void ReferenceCrystalMatchingLog::dump(const CrystalContextResult& crystalContext){
    const AssignmentStats& stats = crystalContext.assignmentStats;
    const double cutoff = crystalContext.selectedCutoff;
    if(crystalContext.tightestCutoffFallback){
        printReciprocityViolations(stats);
        spdlog::warn("No cutoff kept snap p90 below rho {:.3f} A; using the tightest, {:.2f} A",
                     crystalContext.ambiguityRho, cutoff);
    }
    // The selected cutoff's time spans every cutoff tried, so it comes from the selection's own clock; this only
    // ends the phase.
    phaseMilliseconds();
    if(!crystalContext.message.empty()){
        spdlog::error("Bond cutoff selection failed: {}", crystalContext.message);
        return;
    }
    spdlog::info("Selected bond cutoff {:.2f} A in {}: snap p90 {:.3f} A, p99 {:.3f} A (rho {:.3f} A), contract {}",
                 cutoff, _cutoffSelection.elapsed(), crystalContext.snapResidualP90, crystalContext.snapResidualP99,
                 crystalContext.ambiguityRho, crystalContext.contractValid ? "VALID" : "INVALID");

    const std::vector<double>& perAtomResidual = crystalContext.perAtomResidual;
    std::size_t assignedAtoms = 0;
    double residualSum = 0.0;
    double residualMax = 0.0;
    for(double residual : perAtomResidual){
        if(residual < 0.0){
            continue;
        }
        ++assignedAtoms;
        residualSum += residual;
        residualMax = std::max(residualMax, residual);
    }
    const std::size_t atomCount = perAtomResidual.size();
    spdlog::info("  Assigned {} of {} atoms to reference sites ({:.1f}%); rcm_residual mean {:.3f} A, max {:.3f} A",
                 assignedAtoms, atomCount, percentOf(assignedAtoms, atomCount),
                 assignedAtoms > 0 ? residualSum / static_cast<double>(assignedAtoms) : 0.0, residualMax);
    spdlog::info("  Neighbor graph: {} edges, coordination {}-{}, {:.1f}% with ideal lattice vectors",
                 stats.edgesTotal, stats.minNeighborCount, stats.maxNeighborCount,
                 percentOf(stats.edgesWithIdealVector, stats.edgesTotal));
    if(stats.bulkLoopsSampled > 0){
        spdlog::info("  Bulk triangle closure: median {:.2e} A over {} loops",
                     stats.bulkLoopMedianResidual, stats.bulkLoopsSampled);
    }
    if(assignedAtoms < atomCount){
        spdlog::warn("{} atoms left without a reference site ({} with a species absent from the reference, "
                     "{} with no neighbor within {:.2f} A)",
                     atomCount - assignedAtoms, stats.atomsBySpeciesUnassigned, stats.atomsWithZeroNeighbors, cutoff);
    }
}

void ReferenceCrystalMatchingLog::dump(const ReferenceCrystalMatchingService::ContractReport& contract){
    if(!contract.path.empty()){
        if(contract.written){
            spdlog::info("Contract validity report written to {}", contract.path);
        }else{
            spdlog::warn("Could not write the contract validity report to {}", contract.path);
        }
    }
    if(contract.valid){
        spdlog::info("Feed OpenDXA with --reference_topology {}", contract.topologyName);
    }
}

}
