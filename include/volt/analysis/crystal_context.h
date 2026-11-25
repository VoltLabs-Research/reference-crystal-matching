#pragma once

#include <volt/analysis/structure_analysis.h>
#include <volt/analysis/structure_analysis_context.h>
#include <volt/core/lammps_parser.h>
#include <volt/core/particle_property.h>
#include <volt/registration/site_assignment.h>

#include <algorithm>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace Volt{

class ReferenceCrystalMatchingLog;

struct CrystalContextParams{
    std::string referenceFile;
    int anchorSpecies = 1;
    double bondCutoff = 0.0;
    std::string topologyName = "crystal";
};

struct CrystalContextResult{
    bool ok = false;
    std::string message;

    int anchorAtoms = 0;
    double grainSnapResidual = 0;
    double selectedCutoff = 0;

    int grainCount = 1;

    double snapResidualP99 = 0;
    double snapResidualP90 = 0;
    double ambiguityRho = 0;
    bool contractValid = false;
    // No candidate cutoff kept the snap p90 below rho, so the tightest one was used.
    bool tightestCutoffFallback = false;
    AssignmentStats assignmentStats;

    std::shared_ptr<ParticleProperty> structureTypesStorage;

    std::vector<double> perAtomResidual;
};

namespace CrystalContextDetail{

// Half the shortest reference cell edge: a snap residual beyond it no longer tells two sites apart.
inline double ambiguityRho(double cellLengthA, double cellLengthB, double cellLengthC){
    return 0.5 * std::min({cellLengthA, cellLengthB, cellLengthC});
}

// Hands the analysis the neighbor graph, lattice vectors and grains of the chosen assignment, one cluster per grain.
void emitStructureContext(
    StructureContext& context,
    StructureAnalysis& analysis,
    AllAtomNeighbors& neighbors,
    std::vector<Vector3>& overrides,
    std::size_t atomCount,
    double chosenCutoff,
    const std::vector<int>& grainOfAtom,
    const std::vector<Matrix3>& grainRotations,
    const std::string& topologyName,
    CrystalContextResult& result);

}

CrystalContextResult buildCrystalContext(
    const LammpsParser::Frame& frame,
    const std::vector<int>& atomSpecies,
    StructureContext& context,
    StructureAnalysis& analysis,
    const CrystalContextParams& params,
    ReferenceCrystalMatchingLog& log
);

}
