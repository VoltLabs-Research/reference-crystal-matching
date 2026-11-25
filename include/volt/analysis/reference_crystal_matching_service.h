#pragma once

#include <volt/core/volt.h>
#include <nlohmann/json.hpp>
#include <volt/core/lammps_parser.h>
#include <string>
#include <vector>

namespace Volt{

using json = nlohmann::json;

class ReferenceCrystalMatchingService{
public:
    // A reference topology of the lattice directory: the perfect crystal it names and the species anchoring it.
    struct ExoticReference{
        std::string name;
        std::string referenceCrystal;
        int anchorSpecies = 0;
        double fullCrystalCutoff = 0.0;
        bool found = false;
    };

    // Whether the reference registers the sample unambiguously, and where the report saying so went (nowhere
    // without an output base). A run whose contract is not valid fails.
    struct ContractReport{
        std::string path;
        bool written = false;
        bool valid = false;
        std::string topologyName;
    };

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
    // The per-atom species from the requested dump column, falling back to the atom types.
    static bool resolveAtomSpecies(const LammpsParser::Frame& frame,
                                   const std::string& speciesColumn,
                                   std::vector<int>& out,
                                   std::string& error);
    static ExoticReference loadExoticReference(const std::string& latticeDir,
                                               const std::string& name,
                                               std::string& error);

    std::string _referenceCrystal;
    int _anchorSpecies = 1;
    double _fullCrystalCutoff = 0.0;
    std::string _topologyName = "crystal";
    std::string _speciesColumn = "species";
    std::string _referenceTopology;
    std::string _latticeDir;
};

}
