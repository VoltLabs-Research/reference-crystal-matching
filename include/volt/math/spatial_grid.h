#pragma once

#include <volt/math/point3.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <unordered_map>
#include <vector>

namespace Volt{

class SpatialGrid{
public:
    SpatialGrid(const Point3* points, std::size_t count, double cellSize)
        : points_(points), cellSize_(cellSize){
        boundsLow_ = points[0];
        boundsHigh_ = points[0];
        const int n = static_cast<int>(count);
        for(int i = 0; i < n; ++i){
            for(int axis = 0; axis < 3; ++axis){
                boundsLow_[axis]  = std::min(boundsLow_[axis],  points[i][axis]);
                boundsHigh_[axis] = std::max(boundsHigh_[axis], points[i][axis]);
            }
        }
        countX_ = std::max(1, static_cast<int>((boundsHigh_.x() - boundsLow_.x()) / cellSize_) + 1);
        countY_ = std::max(1, static_cast<int>((boundsHigh_.y() - boundsLow_.y()) / cellSize_) + 1);
        countZ_ = std::max(1, static_cast<int>((boundsHigh_.z() - boundsLow_.z()) / cellSize_) + 1);
        for(int i = 0; i < n; ++i){
            const std::array<int, 3> c = cellOf(points[i]);
            grid_[key(c[0], c[1], c[2])].push_back(i);
        }
    }

    const Point3& boundsLow() const{ return boundsLow_; }
    const Point3& boundsHigh() const{ return boundsHigh_; }

    template<typename Fn>
    void forEachNeighbor(int queryIndex, Fn&& fn) const{
        const std::array<int, 3> c = cellOf(points_[queryIndex]);
        for(int dx = -1; dx <= 1; ++dx){
            for(int dy = -1; dy <= 1; ++dy){
                for(int dz = -1; dz <= 1; ++dz){
                    const int cx = c[0] + dx, cy = c[1] + dy, cz = c[2] + dz;
                    if(cx < 0 || cy < 0 || cz < 0 || cx >= countX_ || cy >= countY_ || cz >= countZ_){
                        continue;
                    }
                    const auto bucket = grid_.find(key(cx, cy, cz));
                    if(bucket == grid_.end()){
                        continue;
                    }
                    for(int otherIndex : bucket->second){
                        if(otherIndex != queryIndex){
                            fn(otherIndex);
                        }
                    }
                }
            }
        }
    }

private:
    std::array<int, 3> cellOf(const Point3& p) const{
        return {std::min(countX_ - 1, std::max(0, static_cast<int>((p.x() - boundsLow_.x()) / cellSize_))),
                std::min(countY_ - 1, std::max(0, static_cast<int>((p.y() - boundsLow_.y()) / cellSize_))),
                std::min(countZ_ - 1, std::max(0, static_cast<int>((p.z() - boundsLow_.z()) / cellSize_)))};
    }
    long long key(int x, int y, int z) const{
        return (static_cast<long long>(x) * countY_ + y) * countZ_ + z;
    }

    const Point3* points_;
    double cellSize_;
    Point3 boundsLow_, boundsHigh_;
    int countX_ = 1, countY_ = 1, countZ_ = 1;
    std::unordered_map<long long, std::vector<int>> grid_;
};

}
