#include "grain_frame.h"
#include "statistics.h"

#include <ptm_polar.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace Volt{

namespace{

// An orthonormal frame built from two (non-parallel) vectors by Gram-Schmidt:
// column 0 along the first, column 1 in their plane, column 2 their cross.
Matrix3 orthonormalFrame(const Vector3& firstVector, const Vector3& secondVector){
    const Vector3 axis1 = firstVector / firstVector.length();
    Vector3 axis2 = secondVector - axis1 * secondVector.dot(axis1);
    axis2 /= axis2.length();
    const Vector3 axis3 = axis1.cross(axis2);
    return Matrix3(axis1, axis2, axis3);
}

} // namespace

GrainFrame recoverGrainFrame(
    const Point3* positions,
    const int* species,
    std::size_t atomCount,
    const SimulationCell& cell,
    const AnchorReference& anchorReference)
{
    GrainFrame result;
    constexpr int kNearestNeighbors = 16;
    constexpr double bondCutoff = 6.2;

    std::vector<Point3> anchorPositions;
    for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
        if(species[atomIndex] == anchorReference.targetSpecies){
            anchorPositions.push_back(positions[atomIndex]);
        }
    }
    const int anchorCount = static_cast<int>(anchorPositions.size());
    result.anchorAtomCount = anchorCount;
    if(anchorCount < 50 || anchorReference.idealFractional.empty()){
        result.message = "too few anchor atoms or no ideal reference";
        return result;
    }

    const std::vector<Vector3> idealVectors = buildIdealNeighborVectors(anchorReference);
    if(idealVectors.size() < 4){
        result.message = "ideal reference produced too few neighbour vectors";
        return result;
    }

    Point3 boundsLow = anchorPositions[0];
    Point3 boundsHigh = anchorPositions[0];
    for(const Point3& position : anchorPositions){
        for(int axis = 0; axis < 3; ++axis){
            boundsLow[axis] = std::min(boundsLow[axis], position[axis]);
            boundsHigh[axis] = std::max(boundsHigh[axis], position[axis]);
        }
    }

    const double gridCellSize = anchorReference.physicalBondCutoff;
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
    for(int anchorIndex = 0; anchorIndex < anchorCount; ++anchorIndex){
        const std::array<int, 3> gridCell = gridCellOf(anchorPositions[anchorIndex]);
        grid[gridKey(gridCell[0], gridCell[1], gridCell[2])].push_back(anchorIndex);
    }

    std::vector<std::array<int, kNearestNeighbors>> neighborIndices(static_cast<std::size_t>(anchorCount));
    std::vector<std::array<double, kNearestNeighbors>> neighborDistances(static_cast<std::size_t>(anchorCount));
    for(int anchorIndex = 0; anchorIndex < anchorCount; ++anchorIndex){
        const std::array<int, 3> gridCell = gridCellOf(anchorPositions[anchorIndex]);
        std::vector<std::pair<double, int>> candidates;
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
                        if(otherIndex == anchorIndex){
                            continue;
                        }
                        const double distance =
                            cell.wrapVector(anchorPositions[otherIndex] - anchorPositions[anchorIndex]).length();
                        candidates.emplace_back(distance, otherIndex);
                    }
                }
            }
        }
        std::sort(candidates.begin(), candidates.end());
        for(int slot = 0; slot < kNearestNeighbors; ++slot){
            if(slot < static_cast<int>(candidates.size())){
                neighborDistances[static_cast<std::size_t>(anchorIndex)][slot] = candidates[static_cast<std::size_t>(slot)].first;
                neighborIndices[static_cast<std::size_t>(anchorIndex)][slot] = candidates[static_cast<std::size_t>(slot)].second;
            }else{
                neighborDistances[static_cast<std::size_t>(anchorIndex)][slot] = 1e9;
                neighborIndices[static_cast<std::size_t>(anchorIndex)][slot] = -1;
            }
        }
    }

    const AffineTransformation& cellMatrix = cell.matrix();
    const double cellSpanX = cellMatrix.column(0).length();
    const double cellSpanY = cellMatrix.column(1).length();
    const double midpointZ = (boundsLow.z() + boundsHigh.z()) * 0.5;

    std::vector<int> bulkIndices;
    for(int anchorIndex = 0; anchorIndex < anchorCount; ++anchorIndex){
        const Point3& position = anchorPositions[anchorIndex];
        if(std::abs(position.z() - midpointZ) < 0.45 * (boundsHigh.z() - boundsLow.z()) &&
           position.x() > boundsLow.x() + 12 && position.x() < boundsLow.x() + cellSpanX - 12 &&
           position.y() > boundsLow.y() + 12 && position.y() < boundsLow.y() + cellSpanY - 12){
            bulkIndices.push_back(anchorIndex);
        }
    }
    if(bulkIndices.size() < 20){
        bulkIndices.clear();
        for(int anchorIndex = 0; anchorIndex < anchorCount; ++anchorIndex){
            bulkIndices.push_back(anchorIndex);
        }
    }

    std::vector<Vector3> bulkBonds;
    for(int anchorIndex : bulkIndices){
        for(int slot = 1; slot < kNearestNeighbors; ++slot){
            const int neighbor = neighborIndices[static_cast<std::size_t>(anchorIndex)][slot];
            if(neighbor < 0){
                continue;
            }
            const Vector3 bond = cell.wrapVector(anchorPositions[neighbor] - anchorPositions[anchorIndex]);
            if(bond.length() < bondCutoff){
                bulkBonds.push_back(bond);
            }
        }
    }

    const auto refineAndScore = [&](Matrix3 rotation) -> std::pair<Matrix3, double>{
        for(int iteration = 0; iteration < 6; ++iteration){
            double correlation[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
            int matchedCount = 0;
            for(std::size_t bondIndex = 0; bondIndex < bulkBonds.size(); bondIndex += 7){
                const Vector3 measured = bulkBonds[bondIndex];
                const Vector3 inCrystalFrame = rotation.transposed() * measured;
                double bestDistance = 1e9;
                Vector3 bestIdeal = idealVectors[0];
                for(const Vector3& ideal : idealVectors){
                    const double distance = (ideal - inCrystalFrame).length();
                    if(distance < bestDistance){
                        bestDistance = distance;
                        bestIdeal = ideal;
                    }
                }
                if(bestDistance > 1.0){
                    continue;
                }
                correlation[0] += bestIdeal.x() * measured.x();
                correlation[1] += bestIdeal.x() * measured.y();
                correlation[2] += bestIdeal.x() * measured.z();
                correlation[3] += bestIdeal.y() * measured.x();
                correlation[4] += bestIdeal.y() * measured.y();
                correlation[5] += bestIdeal.y() * measured.z();
                correlation[6] += bestIdeal.z() * measured.x();
                correlation[7] += bestIdeal.z() * measured.y();
                correlation[8] += bestIdeal.z() * measured.z();
                ++matchedCount;
            }
            if(matchedCount < 10){
                break;
            }
            double rotationPart[9];
            double stretchPart[9];
            ptm::polar_decomposition_3x3(correlation, true, rotationPart, stretchPart);
            rotation = Matrix3(
                Vector3(rotationPart[0], rotationPart[3], rotationPart[6]),
                Vector3(rotationPart[1], rotationPart[4], rotationPart[7]),
                Vector3(rotationPart[2], rotationPart[5], rotationPart[8]));
        }

        std::vector<double> residuals;
        for(std::size_t bondIndex = 0; bondIndex < bulkBonds.size(); bondIndex += 11){
            const Vector3 inCrystalFrame = rotation.transposed() * bulkBonds[bondIndex];
            double bestDistance = 1e9;
            for(const Vector3& ideal : idealVectors){
                bestDistance = std::min(bestDistance, (ideal - inCrystalFrame).length());
            }
            residuals.push_back(bestDistance);
        }
        return {rotation, medianOf(residuals)};
    };

    double bestResidual = 1e9;
    int seededCount = 0;
    const std::size_t seedStride = std::max<std::size_t>(1, bulkIndices.size() / 6);
    for(std::size_t bulkPosition = 0; bulkPosition < bulkIndices.size() && seededCount < 6; bulkPosition += seedStride){
        const int anchorIndex = bulkIndices[bulkPosition];

        std::vector<Vector3> measuredBonds;
        for(int slot = 1; slot < kNearestNeighbors; ++slot){
            const int neighbor = neighborIndices[static_cast<std::size_t>(anchorIndex)][slot];
            if(neighbor < 0 || neighborDistances[static_cast<std::size_t>(anchorIndex)][slot] > bondCutoff){
                continue;
            }
            measuredBonds.push_back(cell.wrapVector(anchorPositions[neighbor] - anchorPositions[anchorIndex]));
        }
        if(measuredBonds.size() < 6){
            continue;
        }
        ++seededCount;

        const Vector3 measuredFirst = measuredBonds[0];
        Vector3 measuredSecond(0, 0, 0);
        bool haveSecond = false;
        for(std::size_t bondIndex = 1; bondIndex < measuredBonds.size(); ++bondIndex){
            const Vector3 unit = measuredBonds[bondIndex] / measuredBonds[bondIndex].length();
            if(std::abs(unit.dot(measuredFirst / measuredFirst.length())) < 0.85){
                measuredSecond = measuredBonds[bondIndex];
                haveSecond = true;
                break;
            }
        }
        if(!haveSecond){
            continue;
        }

        const Matrix3 measuredFrame = orthonormalFrame(measuredFirst, measuredSecond);
        const double measuredLengthFirst = measuredFirst.length();
        const double measuredLengthSecond = measuredSecond.length();
        const double measuredAngle = std::acos(std::max(-1.0, std::min(1.0,
            measuredFirst.dot(measuredSecond) / (measuredLengthFirst * measuredLengthSecond))));

        for(std::size_t idealFirst = 0; idealFirst < idealVectors.size(); ++idealFirst){
            if(std::abs(idealVectors[idealFirst].length() - measuredLengthFirst) > 0.25){
                continue;
            }
            for(std::size_t idealSecond = 0; idealSecond < idealVectors.size(); ++idealSecond){
                if(idealSecond == idealFirst){
                    continue;
                }
                if(std::abs(idealVectors[idealSecond].length() - measuredLengthSecond) > 0.25){
                    continue;
                }
                const double idealAngle = std::acos(std::max(-1.0, std::min(1.0,
                    idealVectors[idealFirst].dot(idealVectors[idealSecond]) /
                    (idealVectors[idealFirst].length() * idealVectors[idealSecond].length()))));
                if(std::abs(idealAngle - measuredAngle) > 0.12){
                    continue;
                }
                const Matrix3 seedRotation =
                    measuredFrame * orthonormalFrame(idealVectors[idealFirst], idealVectors[idealSecond]).transposed();
                const auto [refinedRotation, residual] = refineAndScore(seedRotation);
                if(residual < bestResidual){
                    bestResidual = residual;
                    result.rotation = refinedRotation;
                }
            }
        }
        if(bestResidual < 0.45){
            break;
        }
    }

    if(bestResidual > 0.9){
        result.message = "grain-basis recovery failed (residual " + std::to_string(bestResidual) + " A)";
        return result;
    }
    result.bulkSnapResidual = bestResidual;
    result.ok = true;
    return result;
}

}
