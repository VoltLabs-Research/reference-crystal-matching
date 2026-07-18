#include <volt/registration/site_assignment.h>
#include <volt/math/statistics.h>

#include <volt/structures/crystal_structure_types.h>

#include <spdlog/spdlog.h>

#include <algorithm>
#include <cmath>
#include <queue>
#include <unordered_map>
#include <vector>

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
        out_.snapResidualP90 = snapResiduals_.empty()
            ? std::numeric_limits<double>::max()
            : percentileOf(snapResiduals_, 0.90);
        return std::move(out_);
    }

private:
    Vector3& overrideAt(std::size_t atomIndex, int slot){
        return out_.overrides[atomIndex * static_cast<std::size_t>(kNeighborSlots) + static_cast<std::size_t>(slot)];
    }

    const Matrix3& grainRotation(std::size_t atomIndex) const{
        return inputs_.grainRotations[static_cast<std::size_t>(inputs_.grainOfAtom[atomIndex])];
    }
    const Matrix3& grainRotationT(std::size_t atomIndex) const{
        return inputs_.grainRotationsTransposed[static_cast<std::size_t>(inputs_.grainOfAtom[atomIndex])];
    }

    bool snapToSite(int species, const Vector3& predicted,
                    Vector3& outFractional, double& outResidual, int& outSiteIndex) const{
        const auto it = sitesBySpecies_.find(species);
        if(it == sitesBySpecies_.end()){
            return false;
        }
        double best = std::numeric_limits<double>::max();
        for(int siteIndex : it->second){
            const Vector3& siteFrac = inputs_.reference.sites[static_cast<std::size_t>(siteIndex)].fractional;
            Vector3 cand;
            for(int axis = 0; axis < 3; ++axis){
                cand[axis] = siteFrac[axis] + std::round(predicted[axis] - siteFrac[axis]);
            }
            const double residual = (inputs_.reference.cellMatrix * (cand - predicted)).length();
            if(residual < best){
                best = residual;
                outFractional = cand;
                outSiteIndex = siteIndex;
            }
        }
        outResidual = best;
        return best < std::numeric_limits<double>::max();
    }

    int seedSiteForAtom(std::size_t atomIndex) const{
        const AllAtomNeighbors& neighbors = out_.neighbors;
        const int start = neighbors.offsets[atomIndex];
        const int end = neighbors.offsets[atomIndex + 1];
        const int species = inputs_.atomSpecies[atomIndex];
        const auto it = sitesBySpecies_.find(species);
        if(it == sitesBySpecies_.end()){
            return -1;
        }
        int bestSite = -1;
        double bestScore = std::numeric_limits<double>::max();
        for(int siteIndex : it->second){
            const BasisSite& site = inputs_.reference.sites[static_cast<std::size_t>(siteIndex)];
            double score = 0;
            for(int slot = start; slot < end; ++slot){
                const Vector3 crystalDelta =
                    grainRotationT(atomIndex) * neighbors.deltas[static_cast<std::size_t>(slot)];
                const int neighborSpecies =
                    inputs_.atomSpecies[static_cast<std::size_t>(neighbors.indices[static_cast<std::size_t>(slot)])];
                double shellBest = cutoff_ * cutoff_;
                for(const auto& [shellSpecies, shellVector] : site.shell){
                    if(shellSpecies == neighborSpecies){
                        shellBest = std::min(shellBest, (shellVector - crystalDelta).squaredLength());
                    }
                }
                score += shellBest;
            }
            if(score < bestScore){
                bestScore = score;
                bestSite = siteIndex;
            }
        }
        return bestSite;
    }

    void assignBasisSites(){
        const AllAtomNeighbors& neighbors = out_.neighbors;
        AssignmentStats& stats = out_.stats;

        sitesBySpecies_.clear();
        for(int siteIndex = 0; siteIndex < static_cast<int>(inputs_.reference.sites.size()); ++siteIndex){
            sitesBySpecies_[inputs_.reference.sites[static_cast<std::size_t>(siteIndex)].species].push_back(siteIndex);
        }
        const Matrix3 cellInverse = inputs_.reference.cellMatrix.inverse();

        siteFractional_.assign(atomCount_, Vector3::Zero());
        out_.perAtomResidual.assign(atomCount_, -1.0);
        out_.basisSiteOfAtom.assign(atomCount_, -1);

        stats.minNeighborCount = std::numeric_limits<int>::max();
        for(std::size_t atomIndex = 0; atomIndex < atomCount_; ++atomIndex){
            const int neighborCount = neighbors.offsets[atomIndex + 1] - neighbors.offsets[atomIndex];
            stats.minNeighborCount = std::min(stats.minNeighborCount, neighborCount);
            stats.maxNeighborCount = std::max(stats.maxNeighborCount, neighborCount);
            stats.edgesTotal += neighborCount;
            if(neighborCount == 0){
                ++stats.atomsWithZeroNeighbors;
            }
        }
        if(stats.minNeighborCount == std::numeric_limits<int>::max()){
            stats.minNeighborCount = 0;
        }

        struct Entry{
            double residual;
            std::size_t atom;
            Vector3 fractional;
            bool operator>(const Entry& other) const{ return residual > other.residual; }
        };
        std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> frontier;

        auto pushNeighbors = [&](std::size_t atomIndex){
            const int start = neighbors.offsets[atomIndex];
            const int end = neighbors.offsets[atomIndex + 1];
            for(int slot = start; slot < end; ++slot){
                const std::size_t neighbor =
                    static_cast<std::size_t>(neighbors.indices[static_cast<std::size_t>(slot)]);
                if(out_.basisSiteOfAtom[neighbor] >= 0){
                    continue;
                }
                const Vector3 crystalDelta =
                    grainRotationT(atomIndex) * neighbors.deltas[static_cast<std::size_t>(slot)];
                const Vector3 predicted = siteFractional_[atomIndex] + cellInverse * crystalDelta;
                Vector3 snapped;
                double residual;
                int siteIndex;
                if(snapToSite(inputs_.atomSpecies[neighbor], predicted, snapped, residual, siteIndex)){
                    frontier.push({residual, neighbor, snapped});
                }
            }
        };

        std::size_t assignedCount = 0;
        for(std::size_t scan = 0; scan < atomCount_; ++scan){
            if(out_.basisSiteOfAtom[scan] >= 0){
                continue;
            }
            if(neighbors.offsets[scan + 1] - neighbors.offsets[scan] == 0){
                continue;
            }
            const int seedSite = seedSiteForAtom(scan);
            if(seedSite < 0){
                ++stats.atomsBySpeciesUnassigned;
                out_.basisSiteOfAtom[scan] = -2;
                continue;
            }
            siteFractional_[scan] = inputs_.reference.sites[static_cast<std::size_t>(seedSite)].fractional;
            out_.basisSiteOfAtom[scan] = seedSite;
            out_.perAtomResidual[scan] = 0.0;
            ++assignedCount;
            pushNeighbors(scan);

            while(!frontier.empty()){
                const Entry entry = frontier.top();
                frontier.pop();
                if(out_.basisSiteOfAtom[entry.atom] >= 0){
                    continue;
                }
                Vector3 snapped;
                double residual;
                int siteIndex;
                if(!snapToSite(inputs_.atomSpecies[entry.atom], entry.fractional, snapped, residual, siteIndex)){
                    continue;
                }
                siteFractional_[entry.atom] = entry.fractional;
                out_.basisSiteOfAtom[entry.atom] = siteIndex;
                out_.perAtomResidual[entry.atom] = entry.residual;
                ++assignedCount;
                pushNeighbors(entry.atom);
            }
        }

        for(std::size_t atomIndex = 0; atomIndex < atomCount_; ++atomIndex){
            if(out_.basisSiteOfAtom[atomIndex] == -2){
                out_.basisSiteOfAtom[atomIndex] = -1;
            }
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
                const Vector3 crystalFrameDelta = grainRotationT(atomIndex) * neighbors.deltas[static_cast<std::size_t>(slot)];

                const Vector3 measuredFractional = cellInverse * crystalFrameDelta;
                const Vector3 basisFractional = neighborFractional - atomFractional;
                Vector3 idealFractional;
                for(int axis = 0; axis < 3; ++axis){
                    const double imageOffset = std::round(measuredFractional[axis] - basisFractional[axis]);
                    idealFractional[axis] = basisFractional[axis] + imageOffset;
                }
                const Vector3 idealCrystalFrame = cellMatrix * idealFractional;

                overrideAt(atomIndex, slot - start) = idealCrystalFrame;
                if(atomInBulk && snapResiduals_.size() < 20000){
                    const Vector3 idealLabFrame = grainRotation(atomIndex) * idealCrystalFrame;
                    snapResiduals_.push_back((neighbors.deltas[static_cast<std::size_t>(slot)] - idealLabFrame).length());
                }
            }
        }
    }

    void enforceReciprocity(){
        const AllAtomNeighbors& neighbors = out_.neighbors;
        AssignmentStats& stats = out_.stats;
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
                if(inputs_.grainOfAtom[atomIndex] != inputs_.grainOfAtom[static_cast<std::size_t>(neighbor)]){
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
                if(inputs_.grainOfAtom[atomIndex] != inputs_.grainOfAtom[static_cast<std::size_t>(neighbor)]){
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
    std::unordered_map<int, std::vector<int>> sitesBySpecies_;
    std::vector<Vector3> siteFractional_;
};

CutoffAssignment assignForCutoff(const SiteMatchInputs& inputs, double cutoff){
    return SiteAssigner(inputs, cutoff).run();
}

}
