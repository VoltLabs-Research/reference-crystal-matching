#include <volt/registration/neighbor_graph.h>
#include <volt/math/spatial_grid.h>

#include <volt/structures/crystal_structure_types.h>

#include <algorithm>
#include <utility>

namespace Volt{

static constexpr int kNeighborSlots = MAX_NEIGHBORS;

AllAtomNeighbors buildAllAtomNeighbors(
    const Point3* positions,
    std::size_t atomCount,
    const SimulationCell& cell,
    double cutoff)
{
    AllAtomNeighbors result;
    const int count = static_cast<int>(atomCount);
    result.offsets.assign(static_cast<std::size_t>(count) + 1, 0);

    const SpatialGrid grid(positions, atomCount, cutoff);

    std::vector<std::vector<std::pair<int, Vector3>>> perAtomNeighbors(static_cast<std::size_t>(count));
    const double cutoffSquared = cutoff * cutoff;
    for(int atomIndex = 0; atomIndex < count; ++atomIndex){
        std::vector<std::pair<double, std::pair<int, Vector3>>> candidates;
        grid.forEachNeighbor(atomIndex, [&](int otherIndex){
            const Vector3 delta = cell.wrapVector(positions[otherIndex] - positions[atomIndex]);
            const double distanceSquared = delta.squaredLength();
            if(distanceSquared < cutoffSquared && distanceSquared > 0.01){
                candidates.emplace_back(distanceSquared, std::make_pair(otherIndex, delta));
            }
        });
        std::sort(candidates.begin(), candidates.end(), [](const auto& left, const auto& right){
            return left.first < right.first;
        });
        const int keepCount = std::min(static_cast<int>(candidates.size()), kNeighborSlots);
        for(int slot = 0; slot < keepCount; ++slot){
            perAtomNeighbors[static_cast<std::size_t>(atomIndex)].push_back(candidates[static_cast<std::size_t>(slot)].second);
        }
    }

    for(int atomIndex = 0; atomIndex < count; ++atomIndex){
        result.offsets[static_cast<std::size_t>(atomIndex) + 1] =
            result.offsets[static_cast<std::size_t>(atomIndex)] +
            static_cast<int>(perAtomNeighbors[static_cast<std::size_t>(atomIndex)].size());
    }
    result.indices.reserve(static_cast<std::size_t>(result.offsets[static_cast<std::size_t>(count)]));
    result.deltas.reserve(result.indices.capacity());
    for(int atomIndex = 0; atomIndex < count; ++atomIndex){
        for(const auto& [neighbor, delta] : perAtomNeighbors[static_cast<std::size_t>(atomIndex)]){
            result.indices.push_back(neighbor);
            result.deltas.push_back(delta);
        }
    }
    return result;
}

}
