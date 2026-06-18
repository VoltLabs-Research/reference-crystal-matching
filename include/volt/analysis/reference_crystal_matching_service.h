#pragma once

#include <volt/core/volt.h>
#include <nlohmann/json.hpp>
#include <volt/core/lammps_parser.h>
#include <string>

namespace Volt{

using json = nlohmann::json;

class ReferenceCrystalMatchingService{
public:
    ReferenceCrystalMatchingService();

    void setReferenceCrystal(std::string path);
    void setAnchorSpecies(int species);
    void setFullCrystalCutoff(double cutoff);
    void setTopologyName(std::string name);
    void setSpeciesColumn(std::string name);
    void setReferenceTopology(std::string name);
    void setLatticeDir(std::string dir);

    json compute(
        const LammpsParser::Frame& frame,
        const std::string& outputBase,
        const std::string& inputDumpPath
    );

private:
    std::string _referenceCrystal;
    int _anchorSpecies = 1;
    double _fullCrystalCutoff = 0.0;
    std::string _topologyName = "crystal";
    std::string _speciesColumn = "species";
    std::string _referenceTopology;
    std::string _latticeDir;
};

}
