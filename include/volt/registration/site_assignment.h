#pragma once

#include <volt/registration/neighbor_graph.h>
#include <volt/registration/reference_lattice.h>

#include <volt/core/lammps_parser.h>
#include <volt/math/matrix3.h>

#include <limits>
#include <vector>

namespace Volt{

struct AssignmentStats{
    int minNeighborCount = 0;
    int maxNeighborCount = 0;
    long atomsWithZeroNeighbors = 0;
    long edgesTotal = 0;
    long edgesWithIdealVector = 0;
    long atomsBySpeciesUnassigned = 0;
    double bulkLoopMedianResidual = -1;
    int bulkLoopsSampled = 0;
};

struct SiteMatchInputs{
    const LammpsParser::Frame& frame;
    const PerfectReference& reference;
    const std::vector<int>& atomSpecies;
    const std::vector<int>& grainOfAtom;
    const std::vector<Matrix3>& grainRotations;
    const std::vector<Matrix3>& grainRotationsTransposed;
};

struct CutoffAssignment{
    AllAtomNeighbors neighbors;
    std::vector<Vector3> overrides;
    std::vector<int> basisSiteOfAtom;
    std::vector<double> perAtomResidual;
    AssignmentStats stats;
    double snapResidualP99 = std::numeric_limits<double>::max();
    double snapResidualP90 = std::numeric_limits<double>::max();
};

CutoffAssignment assignForCutoff(const SiteMatchInputs& inputs, double cutoff);

}
