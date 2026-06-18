#pragma once

#include <volt/registration/reference_lattice.h>

#include <volt/core/simulation_cell.h>
#include <volt/math/matrix3.h>
#include <volt/math/point3.h>

#include <cstddef>
#include <string>

namespace Volt{

struct GrainFrame{
    bool ok = false;
    std::string message;
    Matrix3 rotation = Matrix3::Identity();
    double bulkSnapResidual = 0;
    int anchorAtomCount = 0;
};

GrainFrame recoverGrainFrame(
    const Point3* positions,
    const int* species,
    std::size_t atomCount,
    const SimulationCell& cell,
    const AnchorReference& anchorReference);

}
