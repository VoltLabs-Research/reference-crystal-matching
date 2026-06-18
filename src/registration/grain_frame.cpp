#include <volt/registration/grain_frame.h>
#include <volt/math/spatial_grid.h>
#include <volt/math/statistics.h>

#include <ptm_polar.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <utility>
#include <vector>

namespace Volt{

static Matrix3 orthonormalFrame(const Vector3& firstVector, const Vector3& secondVector){
    const Vector3 axis1 = firstVector / firstVector.length();
    Vector3 axis2 = secondVector - axis1 * secondVector.dot(axis1);
    axis2 /= axis2.length();
    const Vector3 axis3 = axis1.cross(axis2);
    return Matrix3(axis1, axis2, axis3);
}

static double angleBetween(const Vector3& a, const Vector3& b){
    return std::acos(std::max(-1.0, std::min(1.0, a.dot(b) / (a.length() * b.length()))));
}

class GrainFrameSolver{
public:
    GrainFrameSolver(const Point3* positions, const int* species, std::size_t atomCount,
                     const SimulationCell& cell, const AnchorReference& anchorReference)
        : positions_(positions), species_(species), atomCount_(atomCount),
          cell_(cell), anchorReference_(anchorReference){}

    GrainFrame solve(){
        if(!collectAnchors()){
            return result_;
        }
        buildNeighbours();
        gatherBulkBonds();
        seedSearch();
        return result_;
    }

private:
    static constexpr int kNearestNeighbors = 16;
    static constexpr double bondCutoff = 6.2;

    bool collectAnchors(){
        for(std::size_t atomIndex = 0; atomIndex < atomCount_; ++atomIndex){
            if(species_[atomIndex] == anchorReference_.targetSpecies){
                anchorPositions_.push_back(positions_[atomIndex]);
            }
        }
        anchorCount_ = static_cast<int>(anchorPositions_.size());
        result_.anchorAtomCount = anchorCount_;
        if(anchorCount_ < 50 || anchorReference_.idealFractional.empty()){
            result_.message = "too few anchor atoms or no ideal reference";
            return false;
        }

        idealVectors_ = buildIdealNeighborVectors(anchorReference_);
        if(idealVectors_.size() < 4){
            result_.message = "ideal reference produced too few neighbour vectors";
            return false;
        }
        return true;
    }

    void buildNeighbours(){
        const SpatialGrid grid(anchorPositions_.data(), anchorPositions_.size(),
                               anchorReference_.physicalBondCutoff);
        boundsLow_ = grid.boundsLow();
        boundsHigh_ = grid.boundsHigh();

        neighborIndices_.resize(static_cast<std::size_t>(anchorCount_));
        neighborDistances_.resize(static_cast<std::size_t>(anchorCount_));
        for(int anchorIndex = 0; anchorIndex < anchorCount_; ++anchorIndex){
            std::vector<std::pair<double, int>> candidates;
            grid.forEachNeighbor(anchorIndex, [&](int otherIndex){
                const double distance =
                    cell_.wrapVector(anchorPositions_[otherIndex] - anchorPositions_[anchorIndex]).length();
                candidates.emplace_back(distance, otherIndex);
            });
            std::sort(candidates.begin(), candidates.end());
            auto& distRow = neighborDistances_[static_cast<std::size_t>(anchorIndex)];
            auto& idxRow = neighborIndices_[static_cast<std::size_t>(anchorIndex)];
            for(std::size_t slot = 0; slot < kNearestNeighbors; ++slot){
                if(slot < candidates.size()){
                    distRow[slot] = candidates[slot].first;
                    idxRow[slot] = candidates[slot].second;
                }else{
                    distRow[slot] = 1e9;
                    idxRow[slot] = -1;
                }
            }
        }
    }

    void gatherBulkBonds(){
        const AffineTransformation& cellMatrix = cell_.matrix();
        const double cellSpanX = cellMatrix.column(0).length();
        const double cellSpanY = cellMatrix.column(1).length();
        const double midpointZ = (boundsLow_.z() + boundsHigh_.z()) * 0.5;

        for(int anchorIndex = 0; anchorIndex < anchorCount_; ++anchorIndex){
            const Point3& position = anchorPositions_[anchorIndex];
            if(std::abs(position.z() - midpointZ) < 0.45 * (boundsHigh_.z() - boundsLow_.z()) &&
               position.x() > boundsLow_.x() + 12 && position.x() < boundsLow_.x() + cellSpanX - 12 &&
               position.y() > boundsLow_.y() + 12 && position.y() < boundsLow_.y() + cellSpanY - 12){
                bulkIndices_.push_back(anchorIndex);
            }
        }
        if(bulkIndices_.size() < 20){
            bulkIndices_.clear();
            for(int anchorIndex = 0; anchorIndex < anchorCount_; ++anchorIndex){
                bulkIndices_.push_back(anchorIndex);
            }
        }

        for(int anchorIndex : bulkIndices_){
            for(int slot = 1; slot < kNearestNeighbors; ++slot){
                const int neighbor = neighborIndices_[static_cast<std::size_t>(anchorIndex)][slot];
                if(neighbor < 0){
                    continue;
                }
                const Vector3 bond = cell_.wrapVector(anchorPositions_[neighbor] - anchorPositions_[anchorIndex]);
                if(bond.length() < bondCutoff){
                    bulkBonds_.push_back(bond);
                }
            }
        }
    }

    std::pair<Matrix3, double> refineAndScore(Matrix3 rotation) const{
        for(int iteration = 0; iteration < 6; ++iteration){
            double correlation[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
            int matchedCount = 0;
            for(std::size_t bondIndex = 0; bondIndex < bulkBonds_.size(); bondIndex += 7){
                const Vector3 measured = bulkBonds_[bondIndex];
                const Vector3 inCrystalFrame = rotation.transposed() * measured;
                double bestDistance = 1e9;
                Vector3 bestIdeal = idealVectors_[0];
                for(const Vector3& ideal : idealVectors_){
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
        for(std::size_t bondIndex = 0; bondIndex < bulkBonds_.size(); bondIndex += 11){
            const Vector3 inCrystalFrame = rotation.transposed() * bulkBonds_[bondIndex];
            double bestDistance = 1e9;
            for(const Vector3& ideal : idealVectors_){
                bestDistance = std::min(bestDistance, (ideal - inCrystalFrame).length());
            }
            residuals.push_back(bestDistance);
        }
        return {rotation, medianOf(residuals)};
    }

    void seedSearch(){
        double bestResidual = 1e9;
        int seededCount = 0;
        const std::size_t seedStride = std::max<std::size_t>(1, bulkIndices_.size() / 6);
        for(std::size_t bulkPosition = 0; bulkPosition < bulkIndices_.size() && seededCount < 6; bulkPosition += seedStride){
            const int anchorIndex = bulkIndices_[bulkPosition];

            std::vector<Vector3> measuredBonds;
            for(int slot = 1; slot < kNearestNeighbors; ++slot){
                const int neighbor = neighborIndices_[static_cast<std::size_t>(anchorIndex)][slot];
                if(neighbor < 0 || neighborDistances_[static_cast<std::size_t>(anchorIndex)][slot] > bondCutoff){
                    continue;
                }
                measuredBonds.push_back(cell_.wrapVector(anchorPositions_[neighbor] - anchorPositions_[anchorIndex]));
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
            const double measuredAngle = angleBetween(measuredFirst, measuredSecond);

            for(std::size_t idealFirst = 0; idealFirst < idealVectors_.size(); ++idealFirst){
                if(std::abs(idealVectors_[idealFirst].length() - measuredLengthFirst) > 0.25){
                    continue;
                }
                for(std::size_t idealSecond = 0; idealSecond < idealVectors_.size(); ++idealSecond){
                    if(idealSecond == idealFirst){
                        continue;
                    }
                    if(std::abs(idealVectors_[idealSecond].length() - measuredLengthSecond) > 0.25){
                        continue;
                    }
                    const double idealAngle = angleBetween(idealVectors_[idealFirst], idealVectors_[idealSecond]);
                    if(std::abs(idealAngle - measuredAngle) > 0.12){
                        continue;
                    }
                    const Matrix3 seedRotation =
                        measuredFrame * orthonormalFrame(idealVectors_[idealFirst], idealVectors_[idealSecond]).transposed();
                    const auto [refinedRotation, residual] = refineAndScore(seedRotation);
                    if(residual < bestResidual){
                        bestResidual = residual;
                        result_.rotation = refinedRotation;
                    }
                }
            }
            if(bestResidual < 0.45){
                break;
            }
        }

        if(bestResidual > 0.9){
            result_.message = "grain-basis recovery failed (residual " + std::to_string(bestResidual) + " A)";
            return;
        }
        result_.bulkSnapResidual = bestResidual;
        result_.ok = true;
    }

    const Point3* positions_;
    const int* species_;
    std::size_t atomCount_;
    const SimulationCell& cell_;
    const AnchorReference& anchorReference_;

    std::vector<Point3> anchorPositions_;
    int anchorCount_ = 0;
    std::vector<Vector3> idealVectors_;
    Point3 boundsLow_;
    Point3 boundsHigh_;
    std::vector<std::array<int, kNearestNeighbors>> neighborIndices_;
    std::vector<std::array<double, kNearestNeighbors>> neighborDistances_;
    std::vector<int> bulkIndices_;
    std::vector<Vector3> bulkBonds_;
    GrainFrame result_;
};

GrainFrame recoverGrainFrame(
    const Point3* positions,
    const int* species,
    std::size_t atomCount,
    const SimulationCell& cell,
    const AnchorReference& anchorReference)
{
    return GrainFrameSolver(positions, species, atomCount, cell, anchorReference).solve();
}

}
