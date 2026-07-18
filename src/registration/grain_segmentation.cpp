#include <volt/registration/grain_segmentation.h>
#include <volt/registration/grain_frame.h>
#include <volt/registration/neighbor_graph.h>

#include <spdlog/spdlog.h>
#include <tbb/parallel_for.h>

#include <algorithm>
#include <cmath>
#include <queue>
#include <vector>

namespace Volt{

namespace{

bool sameShell(double lenMeasured, double lenIdeal){
    return std::abs(lenMeasured - lenIdeal) <= 0.15 * lenMeasured;
}

double shellResidualAt(const std::vector<Vector3>& measured, const std::vector<Vector3>& ideal,
                       const Matrix3& frame, int& matchedOut){
    const Matrix3 frameT = frame.transposed();
    double sum = 0;
    int matched = 0;
    for(const Vector3& m : measured){
        const Vector3 inCrystal = frameT * m;
        const double lenM = m.length();
        double best = 1e9;
        for(const Vector3& v : ideal){
            if(!sameShell(lenM, v.length())){
                continue;
            }
            best = std::min(best, (v - inCrystal).length());
        }
        if(best < 1.0){
            sum += best;
            ++matched;
        }
    }
    matchedOut = matched;
    return matched ? sum / matched : 1e9;
}

bool fitsFrame(const std::vector<Vector3>& bonds, const std::vector<Vector3>& ideal,
               const Matrix3& frame, double residualTol){
    if(bonds.size() < 3){
        return false;
    }
    int matched = 0;
    const double residual = shellResidualAt(bonds, ideal, frame, matched);
    const int minMatched = std::max<int>(6, static_cast<int>(0.6 * bonds.size()));
    return matched >= minMatched && residual < residualTol;
}

bool sameOrientation(const Matrix3& a, const Matrix3& b, const std::vector<Vector3>& ideal){
    const Matrix3 rel = a.transposed() * b;
    for(const Vector3& v : ideal){
        const Vector3 rv = rel * v;
        double best = 1e9;
        for(const Vector3& w : ideal){
            best = std::min(best, (rv - w).length());
        }
        if(best > 0.25 * v.length()){
            return false;
        }
    }
    return true;
}

std::vector<Vector3> measuredBondsOf(std::size_t anchorIndex, const AllAtomNeighbors& graph){
    std::vector<Vector3> bonds;
    const int start = graph.offsets[anchorIndex];
    const int end = graph.offsets[anchorIndex + 1];
    for(int slot = start; slot < end; ++slot){
        bonds.push_back(graph.deltas[static_cast<std::size_t>(slot)]);
    }
    return bonds;
}

}

GrainSegmentation segmentGrains(
    const Point3* positions,
    const int* species,
    std::size_t atomCount,
    const SimulationCell& cell,
    const AnchorReference& anchorReference,
    const Matrix3& globalRotation)
{
    GrainSegmentation out;
    out.grainOfAtom.assign(atomCount, -1);

    std::vector<int> atomOfAnchor;
    std::vector<Point3> anchorPositions;
    for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
        if(species[atomIndex] == anchorReference.targetSpecies){
            atomOfAnchor.push_back(static_cast<int>(atomIndex));
            anchorPositions.push_back(positions[atomIndex]);
        }
    }
    const std::size_t anchorCount = anchorPositions.size();
    if(anchorCount < 2){
        out.grainRotations.push_back(globalRotation);
        out.grainCount = 1;
        std::fill(out.grainOfAtom.begin(), out.grainOfAtom.end(), 0);
        return out;
    }

    const AllAtomNeighbors anchorGraph = buildAllAtomNeighbors(
        anchorPositions.data(), anchorCount, cell, AnchorReference::physicalBondCutoff);
    const std::vector<Vector3> idealVectors = buildIdealNeighborVectors(anchorReference);

    std::vector<std::vector<Vector3>> bonds(anchorCount);
    tbb::parallel_for(std::size_t(0), anchorCount, [&](std::size_t a){
        bonds[a] = measuredBondsOf(a, anchorGraph);
    });

    constexpr int kMaxGrains = 64;
    constexpr int kMinGrainAnchors = 50;
    constexpr double residualTol = 0.55;
    const int minGrainSize = std::max<int>(kMinGrainAnchors, static_cast<int>(0.05 * anchorCount));
    std::vector<int> grainOfAnchor(anchorCount, -1);
    std::vector<Matrix3> grainFrames;
    Matrix3 candidate = globalRotation;

    for(int g = 0; g < kMaxGrains; ++g){
        const int newGrainId = static_cast<int>(grainFrames.size());
        std::vector<char> claim(anchorCount, 0);
        tbb::parallel_for(std::size_t(0), anchorCount, [&](std::size_t a){
            if(grainOfAnchor[a] < 0 && fitsFrame(bonds[a], idealVectors, candidate, residualTol)){
                claim[a] = 1;
            }
        });
        int assigned = 0;
        for(std::size_t a = 0; a < anchorCount; ++a){
            if(claim[a] && grainOfAnchor[a] < 0){
                grainOfAnchor[a] = newGrainId;
                ++assigned;
            }
        }
        if(assigned == 0){
            break;
        }
        if(g > 0 && assigned < minGrainSize){
            for(std::size_t a = 0; a < anchorCount; ++a){
                if(grainOfAnchor[a] == newGrainId){
                    grainOfAnchor[a] = -1;
                }
            }
            break;
        }
        grainFrames.push_back(candidate);

        std::vector<Point3> rest;
        for(std::size_t a = 0; a < anchorCount; ++a){
            if(grainOfAnchor[a] < 0){
                rest.push_back(anchorPositions[a]);
            }
        }
        if(static_cast<int>(rest.size()) < kMinGrainAnchors){
            break;
        }
        std::vector<int> restSpecies(rest.size(), anchorReference.targetSpecies);
        const GrainFrame nextFrame = recoverGrainFrame(
            rest.data(), restSpecies.data(), rest.size(), cell, anchorReference);
        if(!nextFrame.ok){
            break;
        }
        bool duplicate = false;
        for(const Matrix3& existing : grainFrames){
            if(sameOrientation(existing, nextFrame.rotation, idealVectors)){
                duplicate = true;
                break;
            }
        }
        if(duplicate){
            break;
        }
        candidate = nextFrame.rotation;
    }

    out.grainCount = static_cast<int>(grainFrames.size());
    if(out.grainCount == 0){
        out.grainRotations.push_back(globalRotation);
        out.grainCount = 1;
        std::fill(out.grainOfAtom.begin(), out.grainOfAtom.end(), 0);
        return out;
    }
    out.grainRotations = std::move(grainFrames);

    for(std::size_t a = 0; a < anchorCount; ++a){
        out.grainOfAtom[static_cast<std::size_t>(atomOfAnchor[a])] = grainOfAnchor[a];
    }
    const AllAtomNeighbors fullGraph = buildAllAtomNeighbors(
        positions, atomCount, cell, AnchorReference::physicalBondCutoff);
    {
        std::queue<std::size_t> flood;
        for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
            if(out.grainOfAtom[atomIndex] >= 0){
                flood.push(atomIndex);
            }
        }
        while(!flood.empty()){
            const std::size_t u = flood.front();
            flood.pop();
            const int start = fullGraph.offsets[u];
            const int end = fullGraph.offsets[u + 1];
            for(int slot = start; slot < end; ++slot){
                const std::size_t v = static_cast<std::size_t>(fullGraph.indices[static_cast<std::size_t>(slot)]);
                if(out.grainOfAtom[v] < 0){
                    out.grainOfAtom[v] = out.grainOfAtom[u];
                    flood.push(v);
                }
            }
        }
        for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
            if(out.grainOfAtom[atomIndex] < 0){
                out.grainOfAtom[atomIndex] = 0;
            }
        }
    }

    spdlog::info("CrystalContext: grain segmentation -> {} grain(s)", out.grainCount);
    return out;
}

}
