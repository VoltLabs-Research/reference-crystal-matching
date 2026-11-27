#include <volt/analysis/reference_crystal_matching_service.h>
#include <volt/analysis/crystal_context.h>
#include <volt/analysis/reconstructed_analysis_pipeline.h>
#include <volt/analysis/reference_crystal_matching_log.h>
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
#include <filesystem>
#include <string>
#include <fstream>
#include <utility>
#include <vector>

namespace Volt{

bool ReferenceCrystalMatchingService::resolveAtomSpecies(const LammpsParser::Frame& frame,
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

ReferenceCrystalMatchingService::ExoticReference ReferenceCrystalMatchingService::loadExoticReference(
                                    const std::string& latticeDir,
                                    const std::string& name,
                                    std::string& error){
    ReferenceCrystalMatchingService::ExoticReference exotic;
    exotic.name = name;
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
    ReferenceCrystalMatchingLog log;
    std::string topologyName = _topologyName;
    if(_referenceCrystal.empty() && !_referenceTopology.empty()){
        std::string error;
        const ExoticReference exotic = loadExoticReference(_latticeDir, _referenceTopology, error);
        if(!exotic.found){
            return AnalysisResult::failure("Reference topology: " + error);
        }
        _referenceCrystal = exotic.referenceCrystal;
        _anchorSpecies = exotic.anchorSpecies;
        if(_fullCrystalCutoff <= 0.0){
            _fullCrystalCutoff = exotic.fullCrystalCutoff;
        }
        topologyName = _referenceTopology;
        log.dump(exotic);
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
    log.dump(atomSpecies);

    try{
        StructureAnalysis analysis(context);

        CrystalContextParams params;
        params.referenceFile = _referenceCrystal;
        params.anchorSpecies = _anchorSpecies;
        params.bondCutoff    = _fullCrystalCutoff;
        params.topologyName  = topologyName;

        // Weighed by how long each step took on a 135k-atom forsterite sample.
        ProgressSteps steps({
            {"Building the crystal context", 8.0},
            {"Writing cluster outputs", 2.0},
            {"Writing the atoms table", 0.5},
        });
        CrystalContextResult fc =
            buildCrystalContext(frame, atomSpecies, context, analysis, params, log);
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

        result["selected_cutoff"] = fc.selectedCutoff;
        result["grain_snap_residual"] = fc.grainSnapResidual;
        result["topology_name"] = topologyName;

        ContractReport contract{{}, false, fc.contractValid, topologyName};
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
                contract.path = outputBase + "_contract_valid.json";
                std::ofstream out(contract.path);
                out << contractValid.dump(2);
                contract.written = static_cast<bool>(out);
            }
            result["contract_valid"] = std::move(contractValid);
        }
        log.dump(contract);
        if(!fc.contractValid){
            return AnalysisResult::failure(
                "RCM contract INVALID: snap p90 " + std::to_string(fc.snapResidualP90) +
                " A >= rho " + std::to_string(fc.ambiguityRho) +
                " A. The reference crystal does not register this sample unambiguously "
                "(incoherent/amorphous/wrong reference). No contract emitted. Override "
                "with an explicit --full_crystal_cutoff only if you know the sample is coherent.");
        }

        steps.next();
        const auto clusterOutputs = AnalysisPipelineUtils::appendClusterOutputs(
            frame, outputBase, annotatedDumpPath, context, analysis, result, &frameError
        );
        if(!clusterOutputs){
            return AnalysisResult::failure(frameError);
        }
        log.dump(clusterOutputs);
        steps.next();

        if(!outputBase.empty()){
            const std::vector<double>& residual = fc.perAtomResidual;
            StructureIdentificationExport::AtomColumnWriter writeResidual;
            if(!residual.empty()){
                writeResidual = [&residual](ColumnarRowWriter& writer,
                                            std::size_t atomIndex, int ){
                    const double r = atomIndex < residual.size() ? residual[atomIndex] : -1.0;
                    writer.field("rcm_residual", r);
                };
            }
            log.dump(StructureIdentificationExport::writeAtomsTable(outputBase, frame, analysis, {}, writeResidual));
        }

        return result;
    }catch(const std::exception& error){
        return AnalysisResult::failure(std::string("Reference-crystal matching failed: ") + error.what());
    }
}

}
