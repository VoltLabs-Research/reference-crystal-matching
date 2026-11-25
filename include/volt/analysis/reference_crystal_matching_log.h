#pragma once

#include <volt/analysis/reference_crystal_matching_service.h>
#include <volt/analysis/structure_identification_log.h>
#include <volt/core/progress.h>

#include <cstddef>
#include <string>
#include <vector>

namespace Volt{

struct AssignmentStats;
struct BasisSite;
struct CrystalContextResult;
struct CutoffAssignment;
struct GrainFrame;
struct GrainSegmentation;
struct PerfectReference;

// The execution log of one reference crystal matching run, on top of the lines every structure identification
// prints.
class ReferenceCrystalMatchingLog : public StructureIdentificationLog{
public:
    using StructureIdentificationLog::dump;

    void dump(const ReferenceCrystalMatchingService::ExoticReference& reference);
    void dump(const std::vector<int>& atomSpecies);
    void dump(const PerfectReference& reference);
    void dump(const GrainFrame& grainFrame);
    // The bond cutoff selection starts once the grains are segmented; its line counts from here.
    void dump(const GrainSegmentation& grains);
    void dump(const CutoffAssignment& trial);
    void dump(const CrystalContextResult& crystalContext);
    void dump(const ReferenceCrystalMatchingService::ContractReport& contract);

private:
    static constexpr std::size_t kMaxListedGrains = 20;

    // "1: 32768 (14.3%), 2: 65536 (28.6%)". Only the first few distinct ids get their own entry, so a wrong
    // --species_column (atom ids, say) cannot blow up the line.
    static std::string describeSpeciesCounts(const std::vector<int>& atomSpecies);
    static std::string describeSiteSpecies(const std::vector<BasisSite>& sites);
    // "147456 (64.3%), 81920 (35.7%)": atoms per grain, listing only the first few grains.
    static std::string describeGrainSizes(const std::vector<int>& grainOfAtom, int grainCount);
    static void printReciprocityViolations(const AssignmentStats& stats);

    Stopwatch _cutoffSelection;
};

}
