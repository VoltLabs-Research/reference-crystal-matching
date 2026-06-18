#pragma once

#include <volt/core/simulation_cell.h>
#include <volt/math/point3.h>
#include <volt/math/vector3.h>

#include <cstddef>
#include <vector>

namespace Volt{

struct AllAtomNeighbors{
    std::vector<int> offsets;
    std::vector<int> indices;
    std::vector<Vector3> deltas;
};

AllAtomNeighbors buildAllAtomNeighbors(
    const Point3* positions,
    std::size_t atomCount,
    const SimulationCell& cell,
    double cutoff);

}
