#include <volt/registration/site_assignment.h>
#include <volt/math/statistics.h>

#include <volt/structures/crystal_structure_types.h>

#include <spdlog/spdlog.h>

#include <algorithm>
#include <cmath>
#include <unordered_map>

namespace Volt{

static constexpr int kNeighborSlots = MAX_NEIGHBORS;

struct BulkRegion{
    Point3 boundsLow;
    double spanX = 0, spanY = 0, spanZ = 0;

    static BulkRegion of(const LammpsParser::Frame& frame, std::size_t atomCount){
        BulkRegion region;
        const AffineTransformation& cellMatrix = frame.simulationCell.matrix();
        region.spanX = cellMatrix.column(0).length();
        region.spanY = cellMatrix.column(1).length();
        region.spanZ = cellMatrix.column(2).length();
        region.boundsLow = frame.positions[0];
        for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
            for(int axis = 0; axis < 3; ++axis){
                region.boundsLow[axis] = std::min(region.boundsLow[axis], frame.positions[atomIndex][axis]);
            }
        }
        return region;
    }

    bool contains(const Point3& position) const{
        return !(position.x() < boundsLow.x() + 0.2 * spanX || position.x() > boundsLow.x() + 0.8 * spanX ||
                 position.y() < boundsLow.y() + 0.2 * spanY || position.y() > boundsLow.y() + 0.8 * spanY ||
                 position.z() < boundsLow.z() + 0.2 * spanZ || position.z() > boundsLow.z() + 0.8 * spanZ);
    }
};

static int slotInto(const AllAtomNeighbors& neighbors, std::size_t fromAtom, int targetAtom){
    const int start = neighbors.offsets[fromAtom];
    const int end = neighbors.offsets[fromAtom + 1];
    for(int slot = start; slot < end; ++slot){
        if(neighbors.indices[static_cast<std::size_t>(slot)] == targetAtom){
            return slot - start;
        }
    }
    return -1;
}

class SiteAssigner{
public:
    SiteAssigner(const SiteMatchInputs& inputs, double cutoff)
        : inputs_(inputs),
          cutoff_(cutoff),
          atomCount_(static_cast<std::size_t>(inputs.frame.natoms)){
        out_.neighbors = buildAllAtomNeighbors(
            inputs_.frame.positions.data(), atomCount_, inputs_.frame.simulationCell, cutoff_);
        out_.overrides.assign(atomCount_ * static_cast<std::size_t>(kNeighborSlots), Vector3::Zero());
        bulk_ = BulkRegion::of(inputs_.frame, atomCount_);
    }

    CutoffAssignment run(){
        assignBasisSites();
        computeIdealOverrides();
        enforceReciprocity();
        sampleLoopResiduals();
        out_.snapResidualP99 = snapResiduals_.empty()
            ? std::numeric_limits<double>::max()
            : percentileOf(snapResiduals_, 0.99);
        return std::move(out_);
    }

private:
    Vector3& overrideAt(std::size_t atomIndex, int slot){
        return out_.overrides[atomIndex * static_cast<std::size_t>(kNeighborSlots) + static_cast<std::size_t>(slot)];
    }

    void assignBasisSites(){
        const PerfectReference& reference = inputs_.reference;
        const AllAtomNeighbors& neighbors = out_.neighbors;
        AssignmentStats& stats = out_.stats;

        std::unordered_map<int, std::vector<int>> sitesBySpecies;
        for(int siteIndex = 0; siteIndex < static_cast<int>(reference.sites.size()); ++siteIndex){
            sitesBySpecies[reference.sites[static_cast<std::size_t>(siteIndex)].species].push_back(siteIndex);
        }

        stats.minNeighborCount = std::numeric_limits<int>::max();

        out_.perAtomResidual.assign(atomCount_, -1.0);
        out_.basisSiteOfAtom.assign(atomCount_, -1);
        for(std::size_t atomIndex = 0; atomIndex < atomCount_; ++atomIndex){
            const int start = neighbors.offsets[atomIndex];
            const int end = neighbors.offsets[atomIndex + 1];
            const int neighborCount = end - start;
            stats.minNeighborCount = std::min(stats.minNeighborCount, neighborCount);
            stats.maxNeighborCount = std::max(stats.maxNeighborCount, neighborCount);
            if(neighborCount == 0){
                ++stats.atomsWithZeroNeighbors;
                continue;
            }

            const int centralSpecies = inputs_.atomSpecies[atomIndex];
            std::vector<Vector3> crystalFrameDeltas(static_cast<std::size_t>(neighborCount));
            for(int slot = 0; slot < neighborCount; ++slot){
                crystalFrameDeltas[static_cast<std::size_t>(slot)] =
                    inputs_.grainRotationTransposed * neighbors.deltas[static_cast<std::size_t>(start + slot)];
            }

            int bestSite = -1;
            double bestScore = std::numeric_limits<double>::max();
            const auto speciesSites = sitesBySpecies.find(centralSpecies);
            if(speciesSites != sitesBySpecies.end()){
                for(int siteIndex : speciesSites->second){
                    const BasisSite& site = reference.sites[static_cast<std::size_t>(siteIndex)];
                    double score = 0;
                    for(int slot = 0; slot < neighborCount; ++slot){
                        const int neighborSpecies =
                            inputs_.atomSpecies[static_cast<std::size_t>(neighbors.indices[static_cast<std::size_t>(start + slot)])];
                        double bestShellDistance = std::numeric_limits<double>::max();
                        for(const auto& [shellSpecies, shellVector] : site.shell){
                            if(shellSpecies != neighborSpecies){
                                continue;
                            }
                            bestShellDistance = std::min(bestShellDistance,
                                (shellVector - crystalFrameDeltas[static_cast<std::size_t>(slot)]).squaredLength());
                        }
                        if(bestShellDistance < std::numeric_limits<double>::max()){
                            score += bestShellDistance;
                        }else{
                            score += cutoff_ * cutoff_;
                        }
                    }
                    if(score < bestScore){
                        bestScore = score;
                        bestSite = siteIndex;
                    }
                }
            }

            stats.edgesTotal += neighborCount;
            if(bestSite < 0){
                ++stats.atomsBySpeciesUnassigned;
                continue;
            }
            out_.basisSiteOfAtom[atomIndex] = bestSite;
            out_.perAtomResidual[atomIndex] = std::sqrt(bestScore / static_cast<double>(neighborCount));
        }
        if(stats.minNeighborCount == std::numeric_limits<int>::max()){
            stats.minNeighborCount = 0;
        }
    }

    void computeIdealOverrides(){
        const PerfectReference& reference = inputs_.reference;
        const AllAtomNeighbors& neighbors = out_.neighbors;
        const Matrix3& cellMatrix = reference.cellMatrix;
        const Matrix3 cellInverse = cellMatrix.inverse();

        for(std::size_t atomIndex = 0; atomIndex < atomCount_; ++atomIndex){
            const int basisSite = out_.basisSiteOfAtom[atomIndex];
            if(basisSite < 0){
                continue;
            }
            const Vector3& atomFractional = reference.sites[static_cast<std::size_t>(basisSite)].fractional;
            const int start = neighbors.offsets[atomIndex];
            const int end = neighbors.offsets[atomIndex + 1];
            const bool atomInBulk = bulk_.contains(inputs_.frame.positions[atomIndex]);

            for(int slot = start; slot < end; ++slot){
                const int neighbor = neighbors.indices[static_cast<std::size_t>(slot)];
                const int neighborBasisSite = out_.basisSiteOfAtom[static_cast<std::size_t>(neighbor)];
                if(neighborBasisSite < 0){
                    continue;
                }
                const Vector3& neighborFractional = reference.sites[static_cast<std::size_t>(neighborBasisSite)].fractional;
                const Vector3 crystalFrameDelta = inputs_.grainRotationTransposed * neighbors.deltas[static_cast<std::size_t>(slot)];

                const Vector3 measuredFractional = cellInverse * crystalFrameDelta;
                const Vector3 basisFractional = neighborFractional - atomFractional;
                Vector3 idealFractional;
                for(int axis = 0; axis < 3; ++axis){
                    const double imageOffset = std::round(measuredFractional[axis] - basisFractional[axis]);
                    idealFractional[axis] = basisFractional[axis] + imageOffset;
                }
                const Vector3 idealLabFrame = inputs_.grainRotation * (cellMatrix * idealFractional);

                overrideAt(atomIndex, slot - start) = idealLabFrame;
                if(atomInBulk && snapResiduals_.size() < 20000){
                    snapResiduals_.push_back((neighbors.deltas[static_cast<std::size_t>(slot)] - idealLabFrame).length());
                }
            }
        }
    }

    void enforceReciprocity(){
        const AllAtomNeighbors& neighbors = out_.neighbors;
        AssignmentStats& stats = out_.stats;
        // Loop bound from the neighbour graph, NOT atomCount_ (= frame.natoms).
        // Do not unify them: this phase's self-check must walk the graph as built.
        const std::size_t atomCount = neighbors.offsets.size() - 1;

        for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
            const int start = neighbors.offsets[atomIndex];
            const int end = neighbors.offsets[atomIndex + 1];
            for(int slot = start; slot < end; ++slot){
                const int neighbor = neighbors.indices[static_cast<std::size_t>(slot)];
                if(neighbor <= static_cast<int>(atomIndex)){
                    continue;
                }
                const int slotForward = slot - start;
                const int slotReverse = slotInto(neighbors, static_cast<std::size_t>(neighbor), static_cast<int>(atomIndex));
                if(slotReverse < 0){
                    continue;
                }
                Vector3& forwardVector = overrideAt(atomIndex, slotForward);
                Vector3& reverseVector = overrideAt(static_cast<std::size_t>(neighbor), slotReverse);
                const bool forwardZero = forwardVector.isZero(1e-9);
                const bool reverseZero = reverseVector.isZero(1e-9);
                if(!forwardZero){
                    reverseVector = -forwardVector;
                }else if(!reverseZero){
                    forwardVector = -reverseVector;
                }
            }
        }

        stats.edgesWithIdealVector = 0;
        long reciprocityViolations = 0;
        for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
            const int start = neighbors.offsets[atomIndex];
            const int end = neighbors.offsets[atomIndex + 1];
            for(int slot = start; slot < end; ++slot){
                const Vector3& forwardVector = overrideAt(atomIndex, slot - start);
                if(forwardVector.isZero(1e-9)){
                    continue;
                }
                ++stats.edgesWithIdealVector;
                const int neighbor = neighbors.indices[static_cast<std::size_t>(slot)];
                const int slotReverse = slotInto(neighbors, static_cast<std::size_t>(neighbor), static_cast<int>(atomIndex));
                if(slotReverse < 0){
                    continue;
                }
                const Vector3& reverseVector = overrideAt(static_cast<std::size_t>(neighbor), slotReverse);
                if(!reverseVector.isZero(1e-9) && !(forwardVector + reverseVector).isZero(1e-6)){
                    ++reciprocityViolations;
                }
            }
        }
        if(reciprocityViolations > 0){
            spdlog::error("CrystalContext: reciprocity self-check FAILED on {} edges (v_AB != -v_BA); "
                          "Burgers circuits will not close", reciprocityViolations);
        }
    }

    void sampleLoopResiduals(){
        const AllAtomNeighbors& neighbors = out_.neighbors;
        std::vector<double> loopResiduals;
        for(std::size_t atomIndex = 0; atomIndex < atomCount_ && loopResiduals.size() < 2000; ++atomIndex){
            if(!bulk_.contains(inputs_.frame.positions[atomIndex])){
                continue;
            }
            const int start = neighbors.offsets[atomIndex];
            const int end = neighbors.offsets[atomIndex + 1];
            for(int firstSlot = start; firstSlot < end; ++firstSlot){
                const int firstNeighbor = neighbors.indices[static_cast<std::size_t>(firstSlot)];
                for(int secondSlot = firstSlot + 1; secondSlot < end; ++secondSlot){
                    const int secondNeighbor = neighbors.indices[static_cast<std::size_t>(secondSlot)];
                    const int slotFirstToSecond = slotInto(neighbors, static_cast<std::size_t>(firstNeighbor), secondNeighbor);
                    if(slotFirstToSecond < 0){
                        continue;
                    }
                    const int slotSecondToAtom = slotInto(neighbors, static_cast<std::size_t>(secondNeighbor), static_cast<int>(atomIndex));
                    if(slotSecondToAtom < 0){
                        continue;
                    }
                    const Vector3 edgeAtomToFirst = overrideAt(atomIndex, firstSlot - start);
                    const Vector3 edgeFirstToSecond = overrideAt(static_cast<std::size_t>(firstNeighbor), slotFirstToSecond);
                    const Vector3 edgeSecondToAtom = overrideAt(static_cast<std::size_t>(secondNeighbor), slotSecondToAtom);
                    if(edgeAtomToFirst.isZero(1e-9) || edgeFirstToSecond.isZero(1e-9) || edgeSecondToAtom.isZero(1e-9)){
                        continue;
                    }
                    loopResiduals.push_back((edgeAtomToFirst + edgeFirstToSecond + edgeSecondToAtom).length());
                    if(loopResiduals.size() >= 2000){
                        break;
                    }
                }
                if(loopResiduals.size() >= 2000){
                    break;
                }
            }
        }
        out_.stats.bulkLoopsSampled = static_cast<int>(loopResiduals.size());
        out_.stats.bulkLoopMedianResidual = loopResiduals.empty() ? -1 : medianOf(loopResiduals);
    }

    const SiteMatchInputs& inputs_;
    double cutoff_;
    std::size_t atomCount_;
    BulkRegion bulk_;
    CutoffAssignment out_;
    std::vector<double> snapResiduals_;
};

CutoffAssignment assignForCutoff(const SiteMatchInputs& inputs, double cutoff){
    return SiteAssigner(inputs, cutoff).run();
}

}
