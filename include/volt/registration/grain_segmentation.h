#pragma once

#include <volt/registration/neighbor_graph.h>
#include <volt/registration/reference_lattice.h>

#include <volt/core/simulation_cell.h>
#include <volt/math/matrix3.h>
#include <volt/math/point3.h>

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace Volt{

struct GrainSegmentation{
    std::vector<int> grainOfAtom;
    std::vector<Matrix3> grainRotations;
    int grainCount = 0;

    // How the search for grains among the anchor atoms went: the anchors each grain it found claimed, the ones
    // left unassigned and why it stopped. No grain found means none fit the global frame.
    struct AnchorSearch{
        std::size_t anchors = 0;
        std::vector<std::size_t> grainAnchors;
        std::size_t unassignedAnchors = 0;
        std::string stopReason;
    };
    // None with fewer than two anchors, which make a single grain.
    std::optional<AnchorSearch> anchorSearch;
    // Atoms no grain reaches through bonds, put in the first grain.
    std::size_t unreachedAtoms = 0;
};

namespace GrainSegmentationDetail{

bool sameShell(double lenMeasured, double lenIdeal);
double shellResidualAt(const std::vector<Vector3>& measured, const std::vector<Vector3>& ideal,
                       const Matrix3& frame, int& matchedOut);
bool fitsFrame(const std::vector<Vector3>& bonds, const std::vector<Vector3>& ideal,
               const Matrix3& frame, double residualTol);
bool sameOrientation(const Matrix3& a, const Matrix3& b, const std::vector<Vector3>& ideal);
std::vector<Vector3> measuredBondsOf(std::size_t anchorIndex, const AllAtomNeighbors& graph);

}

GrainSegmentation segmentGrains(
    const Point3* positions,
    const int* species,
    std::size_t atomCount,
    const SimulationCell& cell,
    const AnchorReference& anchorReference,
    const Matrix3& globalRotation);

}
