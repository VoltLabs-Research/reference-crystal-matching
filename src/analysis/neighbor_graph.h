#pragma once

#include <volt/core/simulation_cell.h>
#include <volt/math/point3.h>
#include <volt/math/vector3.h>

#include <cstddef>
#include <vector>

namespace Volt{

// CSR-style all-atom neighbour lists: for atom i, neighbours live at
// indices/deltas[offsets[i] .. offsets[i+1]). Capped per atom at MAX_NEIGHBORS,
// nearest first.
struct AllAtomNeighbors{
    std::vector<int> offsets;
    std::vector<int> indices;
    std::vector<Vector3> deltas;
};

// Build the all-atom neighbour lists within `cutoff` Angstrom using a uniform
// spatial grid; periodic images resolved via cell.wrapVector.
AllAtomNeighbors buildAllAtomNeighbors(
    const Point3* positions,
    std::size_t atomCount,
    const SimulationCell& cell,
    double cutoff);

}
