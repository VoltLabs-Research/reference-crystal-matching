#include "reference_lattice.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <sstream>

namespace Volt{

PerfectReference buildPerfectReference(const std::string& referenceFile, double cutoff){
    PerfectReference reference;

    std::ifstream input(referenceFile);
    if(!input){
        reference.message = "cannot open reference file '" + referenceFile + "'";
        return reference;
    }

    double boxXLow = 0, boxXHigh = 0;
    double boxYLow = 0, boxYHigh = 0;
    double boxZLow = 0, boxZHigh = 0;
    bool haveX = false, haveY = false, haveZ = false;

    std::vector<int> species;
    std::vector<std::array<double, 3>> positions;
    std::string line;
    bool insideAtomsSection = false;

    while(std::getline(input, line)){
        const std::size_t firstNonSpace = line.find_first_not_of(" \t");
        if(firstNonSpace == std::string::npos){
            continue;
        }
        const std::string trimmed = line.substr(firstNonSpace);

        if(!insideAtomsSection){
            if(trimmed.find("xlo") != std::string::npos){
                std::istringstream stream(trimmed);
                stream >> boxXLow >> boxXHigh;
                haveX = true;
                continue;
            }
            if(trimmed.find("ylo") != std::string::npos){
                std::istringstream stream(trimmed);
                stream >> boxYLow >> boxYHigh;
                haveY = true;
                continue;
            }
            if(trimmed.find("zlo") != std::string::npos){
                std::istringstream stream(trimmed);
                stream >> boxZLow >> boxZHigh;
                haveZ = true;
                continue;
            }
            if(trimmed.rfind("Atoms", 0) == 0){
                insideAtomsSection = true;
                std::getline(input, line);
                continue;
            }
            continue;
        }

        std::istringstream stream(trimmed);
        std::vector<double> fields;
        double field = 0;
        while(stream >> field){
            fields.push_back(field);
        }
        if(fields.size() < 5){
            if(!positions.empty()){
                break;
            }
            continue;
        }

        const int type = static_cast<int>(fields[1]);
        const double x = fields[fields.size() - 3];
        const double y = fields[fields.size() - 2];
        const double z = fields[fields.size() - 1];
        species.push_back(type);
        positions.push_back({x, y, z});
    }

    if(!haveX || !haveY || !haveZ){
        reference.message = "reference file missing box bounds";
        return reference;
    }

    const double lengthA = boxXHigh - boxXLow;
    const double lengthB = boxYHigh - boxYLow;
    const double lengthC = boxZHigh - boxZLow;
    if(lengthA < 1e-6 || lengthB < 1e-6 || lengthC < 1e-6 || positions.empty()){
        reference.message = "degenerate reference cell";
        return reference;
    }

    reference.cellLengthA = lengthA;
    reference.cellLengthB = lengthB;
    reference.cellLengthC = lengthC;

    const double cellLengths[3] = {lengthA, lengthB, lengthC};
    const int siteCount = static_cast<int>(positions.size());
    const int imageRange = static_cast<int>(std::ceil(cutoff / std::min({lengthA, lengthB, lengthC}))) + 1;

    reference.sites.resize(static_cast<std::size_t>(siteCount));
    for(int siteIndex = 0; siteIndex < siteCount; ++siteIndex){
        BasisSite& site = reference.sites[static_cast<std::size_t>(siteIndex)];
        site.species = species[static_cast<std::size_t>(siteIndex)];
        site.fractional = Vector3(
            (positions[siteIndex][0] - boxXLow) / lengthA,
            (positions[siteIndex][1] - boxYLow) / lengthB,
            (positions[siteIndex][2] - boxZLow) / lengthC);

        for(int otherIndex = 0; otherIndex < siteCount; ++otherIndex){
            for(int imageX = -imageRange; imageX <= imageRange; ++imageX){
                for(int imageY = -imageRange; imageY <= imageRange; ++imageY){
                    for(int imageZ = -imageRange; imageZ <= imageRange; ++imageZ){
                        if(siteIndex == otherIndex && imageX == 0 && imageY == 0 && imageZ == 0){
                            continue;
                        }
                        const Vector3 delta(
                            positions[otherIndex][0] + imageX * cellLengths[0] - positions[siteIndex][0],
                            positions[otherIndex][1] + imageY * cellLengths[1] - positions[siteIndex][1],
                            positions[otherIndex][2] + imageZ * cellLengths[2] - positions[siteIndex][2]);
                        const double distance = delta.length();
                        if(distance > 0.1 && distance < cutoff){
                            site.shell.emplace_back(species[static_cast<std::size_t>(otherIndex)], delta);
                        }
                    }
                }
            }
        }
    }

    reference.ok = true;
    return reference;
}

std::vector<Vector3> buildIdealNeighborVectors(const AnchorReference& anchorReference){
    constexpr double neighborCutoff = 6.15;
    const double lengthA = anchorReference.cellLengthA;
    const double lengthB = anchorReference.cellLengthB;
    const double lengthC = anchorReference.cellLengthC;

    std::vector<Vector3> replicatedPositions;
    for(const Vector3& fractional : anchorReference.idealFractional){
        for(int imageX = -2; imageX <= 2; ++imageX){
            for(int imageY = -2; imageY <= 2; ++imageY){
                for(int imageZ = -2; imageZ <= 2; ++imageZ){
                    replicatedPositions.emplace_back(
                        (fractional.x() + imageX) * lengthA,
                        (fractional.y() + imageY) * lengthB,
                        (fractional.z() + imageZ) * lengthC);
                }
            }
        }
    }

    std::vector<Vector3> idealVectors;
    for(const Vector3& fractional : anchorReference.idealFractional){
        const Vector3 center(fractional.x() * lengthA, fractional.y() * lengthB, fractional.z() * lengthC);
        for(const Vector3& replicated : replicatedPositions){
            const Vector3 delta = replicated - center;
            const double distance = delta.length();
            if(distance < 0.1 || distance > neighborCutoff){
                continue;
            }
            bool duplicate = false;
            for(const Vector3& existing : idealVectors){
                if((existing - delta).length() < 0.15){
                    duplicate = true;
                    break;
                }
            }
            if(!duplicate){
                idealVectors.push_back(delta);
            }
        }
    }
    return idealVectors;
}

}
