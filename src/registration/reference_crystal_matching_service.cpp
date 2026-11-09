#include <volt/analysis/reference_crystal_matching_service.h>
#include <volt/analysis/crystal_context.h>
#include <volt/analysis/reconstructed_analysis_pipeline.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/analysis/structure_analysis_context.h>
#include <volt/analysis/structure_identification_export.h>
#include <volt/core/analysis_result.h>
#include <volt/core/progress.h>
#include <volt/structures/crystal_structure_types.h>

#include <yaml-cpp/yaml.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <string>
#include <fstream>
#include <utility>
#include <vector>

namespace Volt{

static bool resolveAtomSpecies(const LammpsParser::Frame& frame,
                        const std::string& speciesColumn,
                        std::vector<int>& out,
                        std::string& error){
    const std::size_t n = static_cast<std::size_t>(frame.natoms);

    const LammpsParser::AtomColumn* column = frame.findAtomProperty(speciesColumn);
    const bool requestedColumnFound = column != nullptr;
    if(!column && speciesColumn != "type"){
        column = frame.findAtomProperty("type");
    }
    if(column && column->size() == n){
        if(!requestedColumnFound){
            spdlog::warn("Species column '{}' not found; using column 'type'", speciesColumn);
        }
        out.assign(n, 0);
        for(std::size_t i = 0; i < n; ++i){
            switch(column->dataType){
                case DataType::Int:    out[i] = column->ints[i]; break;
                case DataType::Int64:  out[i] = static_cast<int>(column->int64s[i]); break;
                case DataType::Double: out[i] = static_cast<int>(std::lround(column->doubles[i])); break;
                case DataType::Void:   break;
            }
        }
        return true;
    }

    if(frame.types.size() == n){
        bool uniform = true;
        for(std::size_t i = 1; i < n && uniform; ++i){
            if(frame.types[i] != frame.types[0]){
                uniform = false;
            }
        }
        const bool hasTypeColumn =
            std::find(frame.atomColumnOrder.begin(), frame.atomColumnOrder.end(), "type")
                != frame.atomColumnOrder.end();
        const bool isDataFile = frame.findHeaderProperty("atom types") != nullptr;
        if(!uniform || hasTypeColumn || isDataFile){
            if(speciesColumn != "type"){
                spdlog::warn("Species column '{}' not usable; using atom types", speciesColumn);
            }
            out = frame.types;
            return true;
        }
    }

    error = "could not resolve per-atom species; dump has no usable '" + speciesColumn +
            "' or 'type' column. Pass --species_column <name>.";
    return false;
}

// "1: 32768 (14.3%), 2: 65536 (28.6%)". Only the first few distinct ids get their own entry so a
// wrong --species_column (e.g. ids) cannot blow up the line or the bookkeeping.
static std::string describeSpeciesCounts(const std::vector<int>& atomSpecies){
    constexpr std::size_t kMaxListedSpecies = 8;
    std::vector<std::pair<int, std::size_t>> counts;
    std::size_t unlistedAtoms = 0;
    for(int species : atomSpecies){
        auto entry = std::find_if(counts.begin(), counts.end(), [species](const auto& item){
            return item.first == species;
        });
        if(entry != counts.end()){
            ++entry->second;
        }else if(counts.size() < kMaxListedSpecies){
            counts.emplace_back(species, 1);
        }else{
            ++unlistedAtoms;
        }
    }
    std::sort(counts.begin(), counts.end());

    std::string text;
    char item[64];
    for(const auto& [species, count] : counts){
        std::snprintf(item, sizeof(item), "%s%d: %zu (%.1f%%)", text.empty() ? "" : ", ",
                      species, count, percentOf(count, atomSpecies.size()));
        text += item;
    }
    if(unlistedAtoms > 0){
        std::snprintf(item, sizeof(item), ", other species: %zu (%.1f%%)",
                      unlistedAtoms, percentOf(unlistedAtoms, atomSpecies.size()));
        text += item;
    }
    return text;
}

struct ExoticReference{
    std::string referenceCrystal;
    int anchorSpecies = 0;
    double fullCrystalCutoff = 0.0;
    bool found = false;
};

static ExoticReference loadExoticReference(const std::string& latticeDir,
                                    const std::string& name,
                                    std::string& error){
    ExoticReference exotic;
    if(latticeDir.empty() || name.empty()){
        return exotic;
    }

    namespace fs = std::filesystem;
    fs::path yamlPath;
    for(const char* extension : {".yml", ".yaml"}){
        const fs::path candidate = fs::path(latticeDir) / (name + extension);
        std::error_code errorCode;
        if(fs::exists(candidate, errorCode)){
            yamlPath = candidate;
            break;
        }
    }
    if(yamlPath.empty()){
        error = "reference topology '" + name + "' not found under '" + latticeDir + "'";
        return exotic;
    }

    YAML::Node document;
    try{
        document = YAML::LoadFile(yamlPath.string());
    }catch(const std::exception& e){
        error = std::string("failed to parse '") + yamlPath.string() + "': " + e.what();
        return exotic;
    }

    if(!document || !document.IsMap()){
        error = "reference topology '" + name + "' is not a YAML map";
        return exotic;
    }
    const YAML::Node speciesNode = document["anchor_species"]
        ? document["anchor_species"]
        : document["cation_species"];
    if(!speciesNode){
        error = "reference topology '" + name + "' is missing anchor_species";
        return exotic;
    }

    exotic.anchorSpecies = speciesNode.as<int>();
    if(!document["reference_crystal"]){
        error = "reference topology '" + name + "' has anchor_species but no reference_crystal";
        return exotic;
    }

    fs::path referencePath(document["reference_crystal"].as<std::string>());
    if(referencePath.is_relative()){
        referencePath = yamlPath.parent_path() / referencePath;
    }
    std::error_code errorCode;
    if(!fs::exists(referencePath, errorCode)){
        error = "reference topology '" + name + "': reference_crystal not found at '" +
                referencePath.string() + "'";
        return exotic;
    }
    exotic.referenceCrystal = referencePath.string();
    if(document["full_crystal_cutoff"]){
        exotic.fullCrystalCutoff = document["full_crystal_cutoff"].as<double>();
    }
    exotic.found = true;
    return exotic;
}

ReferenceCrystalMatchingService::ReferenceCrystalMatchingService() = default;

void ReferenceCrystalMatchingService::setReferenceCrystal(std::string path){ _referenceCrystal = std::move(path); }
void ReferenceCrystalMatchingService::setAnchorSpecies(int species){ _anchorSpecies = species; }
void ReferenceCrystalMatchingService::setFullCrystalCutoff(double cutoff){ _fullCrystalCutoff = cutoff; }
void ReferenceCrystalMatchingService::setTopologyName(std::string name){ _topologyName = std::move(name); }
void ReferenceCrystalMatchingService::setSpeciesColumn(std::string name){ _speciesColumn = std::move(name); }
void ReferenceCrystalMatchingService::setReferenceTopology(std::string name){ _referenceTopology = std::move(name); }
void ReferenceCrystalMatchingService::setLatticeDir(std::string dir){ _latticeDir = std::move(dir); }

json ReferenceCrystalMatchingService::compute(
    const LammpsParser::Frame& frame,
    const std::string& outputBase,
    const std::string& inputDumpPath
){
    std::string topologyName = _topologyName;
    if(_referenceCrystal.empty() && !_referenceTopology.empty()){
        std::string error;
        ExoticReference exotic = loadExoticReference(_latticeDir, _referenceTopology, error);
        if(!exotic.found){
            return AnalysisResult::failure("Reference topology: " + error);
        }
        _referenceCrystal = exotic.referenceCrystal;
        _anchorSpecies = exotic.anchorSpecies;
        if(_fullCrystalCutoff <= 0.0){
            _fullCrystalCutoff = exotic.fullCrystalCutoff;
        }
        topologyName = _referenceTopology;
        spdlog::info("Reference topology '{}' resolved to {} (anchor species {})",
                     _referenceTopology,
                     std::filesystem::path(_referenceCrystal).filename().string(),
                     _anchorSpecies);
    }

    if(_referenceCrystal.empty()){
        return AnalysisResult::failure(
            "--reference_crystal (or --reference_topology in --lattice_dir) is required."
        );
    }

    const std::string annotatedDumpPath = outputBase.empty()
        ? inputDumpPath + ".annotated.dump"
        : outputBase + "_annotated.dump";

    std::string frameError;
    auto session = AnalysisPipelineUtils::prepareAnalysisSession(
        frame, LATTICE_OTHER, &frameError
    );
    if(!session){
        return AnalysisResult::failure(frameError);
    }
    AnalysisContext& context = session->context;

    std::vector<int> atomSpecies;
    std::string speciesError;
    if(!resolveAtomSpecies(frame, _speciesColumn, atomSpecies, speciesError)){
        return AnalysisResult::failure("Species resolution: " + speciesError);
    }
    spdlog::info("Atom species: {}", describeSpeciesCounts(atomSpecies));

    try{
        StructureAnalysis analysis(context);

        CrystalContextParams params;
        params.referenceFile = _referenceCrystal;
        params.anchorSpecies = _anchorSpecies;
        params.bondCutoff    = _fullCrystalCutoff;
        params.topologyName  = topologyName;

        // The crystal context, the cluster outputs and the atoms parquet, weighed by how long each
        // took on a 135k-atom forsterite sample.
        ProgressSteps steps({8.0, 2.0, 0.5});
        CrystalContextResult fc =
            buildCrystalContext(frame, atomSpecies, context, analysis, params);
        if(!fc.ok){
            return AnalysisResult::failure("Reference-crystal context: " + fc.message);
        }

        json result = AnalysisResult::success();
        result["main_listing"] = {
            {"total_atoms", frame.natoms},
            {"anchor_atoms", fc.anchorAtoms},
            {"topology_name", topologyName},
            {"selected_cutoff", fc.selectedCutoff},
            {"grain_snap_residual", fc.grainSnapResidual},
        };

        char rescaleArg[96];
        std::snprintf(rescaleArg, sizeof(rescaleArg), "%.10g,%.10g,%.10g",
                      fc.metricRescaleX, fc.metricRescaleY, fc.metricRescaleZ);
        result["metric_rescale"] = {
            {"sx", fc.metricRescaleX}, {"sy", fc.metricRescaleY}, {"sz", fc.metricRescaleZ}
        };
        result["metric_rescale_arg"] = rescaleArg;
        result["selected_cutoff"] = fc.selectedCutoff;
        result["grain_snap_residual"] = fc.grainSnapResidual;
        result["topology_name"] = topologyName;

        {
            json contractValid = {
                {"valid", fc.contractValid},
                {"snap_residual_p90", fc.snapResidualP90},
                {"snap_residual_p99", fc.snapResidualP99},
                {"ambiguity_rho", fc.ambiguityRho},
                {"grain_snap_residual", fc.grainSnapResidual},
                {"grain_count", fc.grainCount},
                {"selected_cutoff", fc.selectedCutoff},
                {"anchor_atoms", fc.anchorAtoms},
                {"topology_name", topologyName},
            };
            if(!outputBase.empty()){
                const std::string contractPath = outputBase + "_contract_valid.json";
                std::ofstream out(contractPath);
                out << contractValid.dump(2);
                if(out){
                    spdlog::info("Contract validity report written to {}", contractPath);
                }else{
                    spdlog::warn("Could not write the contract validity report to {}", contractPath);
                }
            }
            result["contract_valid"] = std::move(contractValid);
        }
        if(!fc.contractValid){
            return AnalysisResult::failure(
                "RCM contract INVALID: snap p90 " + std::to_string(fc.snapResidualP90) +
                " A >= rho " + std::to_string(fc.ambiguityRho) +
                " A. The reference crystal does not register this sample unambiguously "
                "(incoherent/amorphous/wrong reference). No contract emitted. Override "
                "with an explicit --full_crystal_cutoff only if you know the sample is coherent.");
        }

        spdlog::info("Feed OpenDXA with --reference_topology {} --metric_rescale {}",
                     topologyName, rescaleArg);

        steps.next("Writing cluster outputs");
        spdlog::info("Writing annotated dump, cluster graph and neighbor lattice...");
        const Stopwatch exportTimer;
        if(!AnalysisPipelineUtils::appendClusterOutputs(
            frame, outputBase, annotatedDumpPath, context, analysis, result, &frameError
        )){
            return AnalysisResult::failure(frameError);
        }
        spdlog::info("Cluster outputs written in {}", exportTimer.elapsed());
        steps.next("Writing the atoms table");

        if(!outputBase.empty()){
            const std::vector<double>& residual = fc.perAtomResidual;
            StructureIdentificationExport::AtomColumnWriter writeResidual;
            if(!residual.empty()){
                writeResidual = [&residual](ColumnarAtomWriter& writer,
                                            std::size_t atomIndex, int ){
                    const double r = atomIndex < residual.size() ? residual[atomIndex] : -1.0;
                    writer.field("rcm_residual", r);
                };
            }
            const std::string atomsPath = outputBase + "_atoms.parquet";
            StructureIdentificationExport::streamStructureIdentificationToParquet(
                atomsPath,
                frame,
                analysis,
                {},
                writeResidual
            );
            spdlog::info("Exported atoms data to: {}", atomsPath);
        }

        return result;
    }catch(const std::exception& error){
        return AnalysisResult::failure(std::string("Reference-crystal matching failed: ") + error.what());
    }
}

}
