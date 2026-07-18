#include <volt/analysis/crystal_context.h>

#include <volt/registration/grain_frame.h>
#include <volt/registration/grain_segmentation.h>
#include <volt/registration/reference_lattice.h>
#include <volt/registration/site_assignment.h>

#include <volt/structures/cluster_graph.h>
#include <volt/structures/crystal_structure_types.h>

#include <spdlog/spdlog.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace Volt{

static double ambiguityRho(double cellLengthA, double cellLengthB, double cellLengthC){
    return 0.5 * std::min({cellLengthA, cellLengthB, cellLengthC});
}

static void emitStructureContext(
    StructureContext& context,
    StructureAnalysis& analysis,
    AllAtomNeighbors& neighbors,
    std::vector<Vector3>& overrides,
    std::size_t atomCount,
    double chosenCutoff,
    const std::vector<int>& grainOfAtom,
    const std::vector<Matrix3>& grainRotations,
    const std::string& topologyName,
    CrystalContextResult& result)
{
    auto offsetsProperty = std::make_shared<ParticleProperty>(atomCount + 1, DataType::Int, 1, 0, true);
    std::copy(neighbors.offsets.begin(), neighbors.offsets.end(), offsetsProperty->dataInt());
    auto indicesProperty = std::make_shared<ParticleProperty>(neighbors.indices.size(), DataType::Int, 1, 0, true);
    if(!neighbors.indices.empty()){
        std::copy(neighbors.indices.begin(), neighbors.indices.end(), indicesProperty->dataInt());
    }
    auto countsProperty = std::make_shared<ParticleProperty>(atomCount, DataType::Int, 1, 0, true);
    for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
        countsProperty->setInt(static_cast<int>(atomIndex), neighbors.offsets[atomIndex + 1] - neighbors.offsets[atomIndex]);
    }

    context.neighborOffsets = offsetsProperty;
    context.neighborIndices = indicesProperty;
    context.neighborCounts = countsProperty;

    result.structureTypesStorage = std::make_shared<ParticleProperty>(atomCount, DataType::Int, 1, 0, true);
    constexpr int kSyntheticStructureType = 1;
    for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
        result.structureTypesStorage->setInt(static_cast<int>(atomIndex), kSyntheticStructureType);
    }
    context.structureTypes = result.structureTypesStorage.get();

    const int grainCount = static_cast<int>(grainRotations.size());
    context.atomClusters = std::make_shared<ParticleProperty>(atomCount, DataType::Int, 1, 0, true);
    for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
        context.atomClusters->setInt(static_cast<int>(atomIndex), grainOfAtom[atomIndex] + 1);
    }

    context.maximumNeighborDistance = chosenCutoff;

    analysis.setNeighborLatticeVectorOverrides(std::move(overrides), static_cast<std::size_t>(MAX_NEIGHBORS));

    ClusterGraph& clusterGraph = analysis.clusterGraph();
    std::vector<Cluster*> clusters(static_cast<std::size_t>(grainCount), nullptr);
    for(int grain = 0; grain < grainCount; ++grain){
        Cluster* cluster = clusterGraph.createCluster(0, topologyName, grain + 1);
        cluster->orientation = grainRotations[static_cast<std::size_t>(grain)];
        clusterGraph.createSelfTransition(cluster);
        clusters[static_cast<std::size_t>(grain)] = cluster;
    }
    for(int i = 0; i < grainCount; ++i){
        for(int j = i + 1; j < grainCount; ++j){
            const Matrix3 tmIToJ = grainRotations[static_cast<std::size_t>(j)].transposed() *
                                   grainRotations[static_cast<std::size_t>(i)];
            clusterGraph.createClusterTransition(
                clusters[static_cast<std::size_t>(i)], clusters[static_cast<std::size_t>(j)], tmIToJ, 1);
        }
    }
}

CrystalContextResult buildCrystalContext(
    const LammpsParser::Frame& frame,
    const std::vector<int>& atomSpecies,
    StructureContext& context,
    StructureAnalysis& analysis,
    const CrystalContextParams& params)
{
    CrystalContextResult result;
    const std::size_t atomCount = static_cast<std::size_t>(frame.natoms);

    double referenceCellMin = 0;
    {
        const PerfectReference cellOnly = buildPerfectReference(params.referenceFile, 0.1);
        if(!cellOnly.ok){
            result.message = "perfect reference: " + cellOnly.message;
            return result;
        }
        referenceCellMin = std::min({cellOnly.cellLengthA, cellOnly.cellLengthB, cellOnly.cellLengthC});
    }

    const double referenceReach = std::max(params.bondCutoff, 1.8 * referenceCellMin);
    const PerfectReference reference = buildPerfectReference(params.referenceFile, referenceReach);
    if(!reference.ok){
        result.message = "perfect reference: " + reference.message;
        return result;
    }
    spdlog::info("CrystalContext: reference cell=({:.3f},{:.3f},{:.3f}) basis sites={} refReach={:.2f}",
                 reference.cellLengthA, reference.cellLengthB, reference.cellLengthC,
                 reference.sites.size(), referenceReach);

    result.metricRescaleX = 1.0;
    result.metricRescaleY = 1.0;
    result.metricRescaleZ = 1.0;

    AnchorReference anchorReference;
    anchorReference.cellMatrix = reference.cellMatrix;
    anchorReference.targetSpecies = params.anchorSpecies;
    for(const BasisSite& site : reference.sites){
        if(site.species != params.anchorSpecies){
            continue;
        }
        Vector3 fractional = site.fractional;
        for(int axis = 0; axis < 3; ++axis){
            fractional[axis] -= std::floor(fractional[axis]);
        }
        bool duplicate = false;
        for(const Vector3& existing : anchorReference.idealFractional){
            if((existing - fractional).length() < 0.02){
                duplicate = true;
                break;
            }
        }
        if(!duplicate){
            anchorReference.idealFractional.push_back(fractional);
        }
    }

    const GrainFrame grainFrame = recoverGrainFrame(
        frame.positions.data(), atomSpecies.data(), atomCount, frame.simulationCell, anchorReference);
    result.anchorAtoms = grainFrame.anchorAtomCount;
    if(!grainFrame.ok){
        result.message = "grain frame: " + grainFrame.message;
        return result;
    }
    result.grainSnapResidual = grainFrame.bulkSnapResidual;
    const Matrix3 grainRotation = grainFrame.rotation;
    spdlog::info("CrystalContext: grain frame recovered (residual {:.3f} A, {} anchors)",
                 grainFrame.bulkSnapResidual, grainFrame.anchorAtomCount);

    const double rho = ambiguityRho(reference.cellLengthA, reference.cellLengthB, reference.cellLengthC);

    double firstNeighborDist = 0.0;
    for(const BasisSite& site : reference.sites){
        for(const auto& [shellSpecies, shellVector] : site.shell){
            const double d = shellVector.length();
            if(firstNeighborDist <= 0.0 || d < firstNeighborDist){
                firstNeighborDist = d;
            }
        }
    }

    const GrainSegmentation grains = segmentGrains(
        frame.positions.data(), atomSpecies.data(), atomCount, frame.simulationCell,
        anchorReference, grainRotation);
    std::vector<Matrix3> grainRotationsT(grains.grainRotations.size());
    for(std::size_t g = 0; g < grains.grainRotations.size(); ++g){
        grainRotationsT[g] = grains.grainRotations[g].transposed();
    }
    result.grainCount = grains.grainCount;

    const SiteMatchInputs matchInputs{
        frame, reference, atomSpecies, grains.grainOfAtom, grains.grainRotations, grainRotationsT};

    std::vector<double> candidates;
    if(params.bondCutoff > 0.0){
        candidates.push_back(params.bondCutoff);
    }else{
        const double cellMin = std::min({reference.cellLengthA, reference.cellLengthB, reference.cellLengthC});
        const double step = std::max(0.1, cellMin / 24.0);
        const double firstShell = firstNeighborDist > 0.1 ? firstNeighborDist : 0.5 * cellMin;
        const double lowestCutoff = std::max(step, 1.15 * firstShell);
        const double highestCutoff = std::max({0.85 * cellMin, 1.4 * firstShell, lowestCutoff});
        for(double cutoff = highestCutoff; cutoff >= lowestCutoff - 1e-9; cutoff -= step){
            candidates.push_back(cutoff);
        }
    }

    CutoffAssignment assignment;
    double chosenCutoff = candidates.empty() ? params.bondCutoff : candidates.back();
    bool picked = false;
    for(double cutoff : candidates){
        CutoffAssignment trial = assignForCutoff(matchInputs, cutoff);
        spdlog::info("CrystalContext: trial cutoff {:.2f} -> snap p90 {:.3f}A p99 {:.3f}A (rho {:.3f}), coord<= {}, idealEdges {}",
                     cutoff, trial.snapResidualP90, trial.snapResidualP99, rho,
                     trial.stats.maxNeighborCount, trial.stats.edgesWithIdealVector);
        if(trial.snapResidualP90 < rho || params.bondCutoff > 0.0){
            chosenCutoff = cutoff;
            assignment = std::move(trial);
            picked = true;
            break;
        }
        assignment = std::move(trial);
    }
    if(!picked && !candidates.empty()){
        chosenCutoff = candidates.back();
        assignment = assignForCutoff(matchInputs, chosenCutoff);
        spdlog::warn("CrystalContext: no cutoff kept snap p90 < rho {:.3f}A; using tightest {:.2f}A", rho, chosenCutoff);
    }

    if(chosenCutoff <= 0.0 || assignment.neighbors.offsets.size() < atomCount + 1){
        result.selectedCutoff = chosenCutoff;
        result.ambiguityRho = rho;
        result.contractValid = false;
        result.ok = true;
        result.message = "no usable neighbour cutoff for this reference/sample "
                         "(cell too small or no coherent crystal); contract rejected";
        spdlog::error("CrystalContext: {}", result.message);
        return result;
    }

    result.perAtomResidual = std::move(assignment.perAtomResidual);
    result.selectedCutoff = chosenCutoff;
    result.snapResidualP99 = assignment.snapResidualP99;
    result.snapResidualP90 = assignment.snapResidualP90;
    result.ambiguityRho = rho;
    result.contractValid = (assignment.snapResidualP90 < rho) || (params.bondCutoff > 0.0);
    spdlog::info("CrystalContext: SELECTED cutoff {:.2f}A (rho {:.3f}A, snap p90 {:.3f}A p99 {:.3f}A, grains {}, contract {})",
                 chosenCutoff, rho, assignment.snapResidualP90, assignment.snapResidualP99,
                 result.grainCount, result.contractValid ? "VALID" : "INVALID");

    emitStructureContext(context, analysis, assignment.neighbors, assignment.overrides,
                         atomCount, chosenCutoff, grains.grainOfAtom, grains.grainRotations,
                         params.topologyName, result);

    result.ok = true;
    return result;
}

}
