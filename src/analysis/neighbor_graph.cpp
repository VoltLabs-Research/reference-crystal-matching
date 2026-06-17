#include "neighbor_graph.h"

#include <volt/structures/crystal_structure_types.h>

#include <algorithm>
#include <array>
#include <unordered_map>
#include <utility>

namespace Volt{

namespace{
constexpr int kNeighborSlots = MAX_NEIGHBORS;
}

AllAtomNeighbors buildAllAtomNeighbors(
    const Point3* positions,
    std::size_t atomCount,
    const SimulationCell& cell,
    double cutoff)
{
    AllAtomNeighbors result;
    const int count = static_cast<int>(atomCount);
    result.offsets.assign(static_cast<std::size_t>(count) + 1, 0);

    Point3 boundsLow = positions[0];
    Point3 boundsHigh = positions[0];
    for(int atomIndex = 0; atomIndex < count; ++atomIndex){
        for(int axis = 0; axis < 3; ++axis){
            boundsLow[axis] = std::min(boundsLow[axis], positions[atomIndex][axis]);
            boundsHigh[axis] = std::max(boundsHigh[axis], positions[atomIndex][axis]);
        }
    }

    const double gridCellSize = cutoff;
    const int gridCountX = std::max(1, static_cast<int>((boundsHigh.x() - boundsLow.x()) / gridCellSize) + 1);
    const int gridCountY = std::max(1, static_cast<int>((boundsHigh.y() - boundsLow.y()) / gridCellSize) + 1);
    const int gridCountZ = std::max(1, static_cast<int>((boundsHigh.z() - boundsLow.z()) / gridCellSize) + 1);

    const auto gridCellOf = [&](const Point3& position){
        return std::array<int, 3>{
            std::min(gridCountX - 1, std::max(0, static_cast<int>((position.x() - boundsLow.x()) / gridCellSize))),
            std::min(gridCountY - 1, std::max(0, static_cast<int>((position.y() - boundsLow.y()) / gridCellSize))),
            std::min(gridCountZ - 1, std::max(0, static_cast<int>((position.z() - boundsLow.z()) / gridCellSize)))};
    };
    const auto gridKey = [&](int cellX, int cellY, int cellZ){
        return (static_cast<long long>(cellX) * gridCountY + cellY) * gridCountZ + cellZ;
    };

    std::unordered_map<long long, std::vector<int>> grid;
    for(int atomIndex = 0; atomIndex < count; ++atomIndex){
        const std::array<int, 3> gridCell = gridCellOf(positions[atomIndex]);
        grid[gridKey(gridCell[0], gridCell[1], gridCell[2])].push_back(atomIndex);
    }

    std::vector<std::vector<std::pair<int, Vector3>>> perAtomNeighbors(static_cast<std::size_t>(count));
    const double cutoffSquared = cutoff * cutoff;
    for(int atomIndex = 0; atomIndex < count; ++atomIndex){
        const std::array<int, 3> gridCell = gridCellOf(positions[atomIndex]);
        std::vector<std::pair<double, std::pair<int, Vector3>>> candidates;
        for(int offsetX = -1; offsetX <= 1; ++offsetX){
            for(int offsetY = -1; offsetY <= 1; ++offsetY){
                for(int offsetZ = -1; offsetZ <= 1; ++offsetZ){
                    const int cellX = gridCell[0] + offsetX;
                    const int cellY = gridCell[1] + offsetY;
                    const int cellZ = gridCell[2] + offsetZ;
                    if(cellX < 0 || cellY < 0 || cellZ < 0 ||
                       cellX >= gridCountX || cellY >= gridCountY || cellZ >= gridCountZ){
                        continue;
                    }
                    const auto bucket = grid.find(gridKey(cellX, cellY, cellZ));
                    if(bucket == grid.end()){
                        continue;
                    }
                    for(int otherIndex : bucket->second){
                        if(otherIndex == atomIndex){
                            continue;
                        }
                        const Vector3 delta = cell.wrapVector(positions[otherIndex] - positions[atomIndex]);
                        const double distanceSquared = delta.squaredLength();
                        if(distanceSquared < cutoffSquared && distanceSquared > 0.01){
                            candidates.emplace_back(distanceSquared, std::make_pair(otherIndex, delta));
                        }
                    }
                }
            }
        }
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
