#pragma once

#include <volt/registration/reference_lattice.h>

#include <volt/core/simulation_cell.h>
#include <volt/math/matrix3.h>
#include <volt/math/point3.h>

#include <cstddef>
#include <vector>

namespace Volt{

struct GrainSegmentation{
    std::vector<int> grainOfAtom;
    std::vector<Matrix3> grainRotations;
    int grainCount = 0;
};

GrainSegmentation segmentGrains(
    const Point3* positions,
    const int* species,
    std::size_t atomCount,
    const SimulationCell& cell,
    const AnchorReference& anchorReference,
    const Matrix3& globalRotation);

}
