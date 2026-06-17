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
    double tiltXY = 0, tiltXZ = 0, tiltYZ = 0;   // non-orthogonal cell support
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
            if(trimmed.find("xy") != std::string::npos && trimmed.find("xz") != std::string::npos){
                std::istringstream stream(trimmed);
                stream >> tiltXY >> tiltXZ >> tiltYZ;
                continue;
            }
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

    // LAMMPS triclinic convention: cell edge vectors as matrix columns
    //   a = (lengthA, 0, 0), b = (xy, lengthB, 0), c = (xz, yz, lengthC).
    // For an orthogonal cell the tilts are 0 and this is diagonal.
    const Matrix3 cellMatrix(
        Vector3(lengthA, 0.0, 0.0),
        Vector3(tiltXY, lengthB, 0.0),
        Vector3(tiltXZ, tiltYZ, lengthC));
    reference.cellMatrix = cellMatrix;
    const Matrix3 cellInverse = cellMatrix.inverse();

    const int siteCount = static_cast<int>(positions.size());
    const int imageRange = static_cast<int>(std::ceil(cutoff / std::min({lengthA, lengthB, lengthC}))) + 1;

    reference.sites.resize(static_cast<std::size_t>(siteCount));
    for(int siteIndex = 0; siteIndex < siteCount; ++siteIndex){
        BasisSite& site = reference.sites[static_cast<std::size_t>(siteIndex)];
        site.species = species[static_cast<std::size_t>(siteIndex)];
        const Vector3 cart(positions[siteIndex][0] - boxXLow,
                           positions[siteIndex][1] - boxYLow,
                           positions[siteIndex][2] - boxZLow);
        site.fractional = cellInverse * cart;

        for(int otherIndex = 0; otherIndex < siteCount; ++otherIndex){
            const Vector3 base(positions[otherIndex][0] - positions[siteIndex][0],
                               positions[otherIndex][1] - positions[siteIndex][1],
                               positions[otherIndex][2] - positions[siteIndex][2]);
            for(int imageX = -imageRange; imageX <= imageRange; ++imageX){
                for(int imageY = -imageRange; imageY <= imageRange; ++imageY){
                    for(int imageZ = -imageRange; imageZ <= imageRange; ++imageZ){
                        if(siteIndex == otherIndex && imageX == 0 && imageY == 0 && imageZ == 0){
                            continue;
                        }
                        const Vector3 delta = base + cellMatrix *
                            Vector3(static_cast<double>(imageX),
                                    static_cast<double>(imageY),
                                    static_cast<double>(imageZ));
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
    const Matrix3& cellMatrix = anchorReference.cellMatrix;

    // Cartesian position of a fractional coord = cellMatrix * fractional, which
    // reduces to per-axis scaling for an orthogonal cell and handles tilt
    // (HCP/A7 primitive cells) for non-orthogonal ones.
    std::vector<Vector3> replicatedPositions;
    for(const Vector3& fractional : anchorReference.idealFractional){
        for(int imageX = -2; imageX <= 2; ++imageX){
            for(int imageY = -2; imageY <= 2; ++imageY){
                for(int imageZ = -2; imageZ <= 2; ++imageZ){
                    replicatedPositions.push_back(cellMatrix * Vector3(
                        fractional.x() + imageX,
                        fractional.y() + imageY,
                        fractional.z() + imageZ));
                }
            }
        }
    }

    std::vector<Vector3> idealVectors;
    for(const Vector3& fractional : anchorReference.idealFractional){
        const Vector3 center = cellMatrix * fractional;
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
