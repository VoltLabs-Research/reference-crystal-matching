#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace Volt{

inline double medianOf(std::vector<double> values){
    if(values.empty()){
        return 0.0;
    }
    std::sort(values.begin(), values.end());
    return values[values.size() / 2];
}

inline double percentileOf(std::vector<double> values, double fraction){
    if(values.empty()){
        return 0.0;
    }
    std::sort(values.begin(), values.end());
    const double position = fraction * (static_cast<double>(values.size()) - 1.0);
    const std::size_t lowerIndex = static_cast<std::size_t>(std::floor(position));
    const std::size_t upperIndex = static_cast<std::size_t>(std::ceil(position));
    if(lowerIndex == upperIndex){
        return values[lowerIndex];
    }
    const double weight = position - static_cast<double>(lowerIndex);
    return values[lowerIndex] * (1.0 - weight) + values[upperIndex] * weight;
}

}
