#pragma once

#include "reference_lattice.h"

#include <volt/core/simulation_cell.h>
#include <volt/math/matrix3.h>
#include <volt/math/point3.h>

#include <cstddef>
#include <string>

namespace Volt{

// Result of recovering the grain's crystallographic orientation: the rotation
// that maps reference (crystal) frame -> lab frame, plus the snap residual and
// the anchor-atom count that fed the fit.
struct GrainFrame{
    bool ok = false;
    std::string message;
    Matrix3 rotation = Matrix3::Identity();
    double bulkSnapResidual = 0;
    int anchorAtomCount = 0;
};

// Recover the grain orientation by matching the anchor sublattice's measured
// bonds against the ideal anchor neighbour vectors: seed from a few bulk atoms,
// refine each seed by polar decomposition (ICP), keep the lowest-residual fit.
GrainFrame recoverGrainFrame(
    const Point3* positions,
    const int* species,
    std::size_t atomCount,
    const SimulationCell& cell,
    const AnchorReference& anchorReference);

}
