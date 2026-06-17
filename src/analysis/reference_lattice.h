#pragma once

#include <volt/math/matrix3.h>
#include <volt/math/vector3.h>

#include <string>
#include <utility>
#include <vector>

namespace Volt{

// One site of the perfect reference crystal's basis, with its in-cell shell of
// (species, vector) neighbours used for per-atom site assignment.
struct BasisSite{
    int species = 0;
    Vector3 fractional = Vector3::Zero();
    std::vector<std::pair<int, Vector3>> shell;
};

// The perfect reference crystal parsed from a LAMMPS data file: the full cell
// matrix (columns = cell edge vectors a,b,c, including tilt for non-orthogonal
// cells) + every basis site with its precomputed neighbour shell. cellLength*
// are the edge magnitudes, kept for logging / rho heuristics.
struct PerfectReference{
    bool ok = false;
    std::string message;
    Matrix3 cellMatrix = Matrix3(Matrix3::Identity{});
    double cellLengthA = 0;
    double cellLengthB = 0;
    double cellLengthC = 0;
    std::vector<BasisSite> sites;
};

// The anchor sublattice: the ideal fractional positions of the single species
// whose dense, near-isotropic sublattice fixes the grain orientation. (Any
// species works; "anchor" is the role, not a chemistry claim.)
struct AnchorReference{
    Matrix3 cellMatrix = Matrix3(Matrix3::Identity{});
    double cellLengthA = 0;
    double cellLengthB = 0;
    double cellLengthC = 0;
    std::vector<Vector3> idealFractional;
    int targetSpecies = 0;
    static constexpr double physicalBondCutoff = 5.7;
};

// Parse a LAMMPS reference data file into cell lengths + basis sites, computing
// each site's neighbour shell out to `cutoff` Angstrom.
PerfectReference buildPerfectReference(const std::string& referenceFile, double cutoff);

// The set of distinct ideal anchor->anchor neighbour vectors (lab-frame,
// deduplicated) that the grain-frame solver snaps measured bonds against.
std::vector<Vector3> buildIdealNeighborVectors(const AnchorReference& anchorReference);

}
