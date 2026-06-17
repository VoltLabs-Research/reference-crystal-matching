#pragma once

#include <volt/core/volt.h>
#include <nlohmann/json.hpp>
#include <volt/core/lammps_parser.h>
#include <string>

namespace Volt{

using json = nlohmann::json;

// ReferenceCrystalMatching (RCM) producer service.
//
// Registers an observed atomic configuration against a perfect reference
// crystal of arbitrary composition/symmetry (ICP / polar decomposition), then
// assigns every atom its ideal neighbor vectors IN THE REFERENCE LATTICE. Because
// all ideal vectors derive from one shared reference, the resulting elastic
// mapping is reciprocal (v_AB = -v_BA) and closes Burgers circuits — which a
// per-atom local classifier (PTM/CNA) cannot guarantee for multi-species,
// low-symmetry crystals. Emits the structure-identification contract (annotated
// dump + clusters table + cluster transitions table + neighbor-lattice sidecar)
// for OpenDXA to consume.
class ReferenceCrystalMatchingService{
public:
    ReferenceCrystalMatchingService();

    // Perfect reference crystal (.lmp). Required unless provided via
    // --reference_topology + --lattice_dir.
    void setReferenceCrystal(std::string path);
    // Species id (in the reference + dump) whose dense sublattice anchors the
    // grain orientation.
    void setAnchorSpecies(int species);
    // Bond cutoff in Angstrom for the reference shell; 0 = auto-select.
    void setFullCrystalCutoff(double cutoff);
    // Topology name written into the cluster table (consumed by OpenDXA's
    // --reference_topology).
    void setTopologyName(std::string name);
    // Dump column carrying the atomic species. LAMMPS only maps a column named
    // "type" to frame.types; dumps using "species" need this override.
    void setSpeciesColumn(std::string name);
    // Optional convenience: an "exotic" YAML (anchor_species + reference_crystal
    // [+ full_crystal_cutoff]) resolved by basename from --lattice_dir.
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
