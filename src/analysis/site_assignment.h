#pragma once

#include "neighbor_graph.h"
#include "reference_lattice.h"

#include <volt/core/lammps_parser.h>
#include <volt/math/matrix3.h>

#include <limits>
#include <vector>

namespace Volt{

// Per-cutoff quality counters. Previously these lived in FullCrystalContextResult
// and were reused as a scratchpad across trial cutoffs; pulling them out lets the
// result struct mean only "the final answer".
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

// Read-only inputs shared by every assignment phase. A call-local bundle: it
// must not outlive the frame/reference it borrows.
struct SiteMatchInputs{
    const LammpsParser::Frame& frame;
    const PerfectReference& reference;
    const std::vector<int>& atomSpecies;
    Matrix3 grainRotation;
    Matrix3 grainRotationTransposed;
};

// Everything one cutoff produces: the neighbour graph, the per-edge reciprocal
// ideal-vector overrides, the per-atom basis-site assignment, quality stats, and
// the p99 snap residual used to accept/reject the cutoff.
struct CutoffAssignment{
    AllAtomNeighbors neighbors;
    std::vector<Vector3> overrides;     // atomCount * MAX_NEIGHBORS, row-major
    std::vector<int> basisSiteOfAtom;   // atomCount; -1 = unassigned
    // Per-atom RMS deviation (Angstrom) of the atom's neighbour shell from the
    // best-matching reference site. Low = bulk crystal; high = surface or
    // dislocation core. This is RCM's continuous "defectness" signal (the
    // analogue of PTM/CNA's OTHER classification). -1 if unassigned.
    std::vector<double> perAtomResidual;
    AssignmentStats stats;
    double snapResidualP99 = std::numeric_limits<double>::max();
};

// Run the full per-cutoff pipeline: build neighbours, assign each atom to its
// reference basis site, derive reciprocal ideal neighbour vectors, and sample
// bulk loop residuals. The reciprocity invariant (v_AB = -v_BA) is what makes
// Burgers circuits close, so it is enforced and self-checked here.
CutoffAssignment assignForCutoff(const SiteMatchInputs& inputs, double cutoff);

}
