#pragma once

#include <volt/math/matrix3.h>
#include <volt/math/vector3.h>

#include <string>
#include <utility>
#include <vector>

namespace Volt{

struct BasisSite{
    int species = 0;
    Vector3 fractional = Vector3::Zero();
    std::vector<std::pair<int, Vector3>> shell;
};

struct PerfectReference{
    bool ok = false;
    std::string message;
    Matrix3 cellMatrix = Matrix3(Matrix3::Identity{});
    double cellLengthA = 0;
    double cellLengthB = 0;
    double cellLengthC = 0;
    std::vector<BasisSite> sites;
};

struct AnchorReference{
    Matrix3 cellMatrix = Matrix3(Matrix3::Identity{});
    std::vector<Vector3> idealFractional;
    int targetSpecies = 0;
    static constexpr double physicalBondCutoff = 5.7;
};

PerfectReference buildPerfectReference(const std::string& referenceFile, double cutoff);

std::vector<Vector3> buildIdealNeighborVectors(const AnchorReference& anchorReference);

}
