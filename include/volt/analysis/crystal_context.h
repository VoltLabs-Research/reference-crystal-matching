#pragma once

#include <volt/analysis/structure_analysis.h>
#include <volt/analysis/structure_analysis_context.h>
#include <volt/core/lammps_parser.h>
#include <volt/core/particle_property.h>

#include <memory>
#include <string>
#include <vector>

namespace Volt{

struct CrystalContextParams{
    std::string referenceFile;
    int anchorSpecies = 1;
    double bondCutoff = 0.0;
    std::string topologyName = "crystal";
};

struct CrystalContextResult{
    bool ok = false;
    std::string message;

    int anchorAtoms = 0;
    double grainSnapResidual = 0;
    double selectedCutoff = 0;

    double metricRescaleX = 1.0, metricRescaleY = 1.0, metricRescaleZ = 1.0;

    std::shared_ptr<ParticleProperty> structureTypesStorage;

    std::vector<double> perAtomResidual;
};

CrystalContextResult buildCrystalContext(
    const LammpsParser::Frame& frame,
    const std::vector<int>& atomSpecies,
    StructureContext& context,
    StructureAnalysis& analysis,
    const CrystalContextParams& params
);

}
