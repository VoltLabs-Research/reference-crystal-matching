#include <volt/plugin/plugin_entry.h>
#include <volt/analysis/reference_crystal_matching_service.h>

using namespace Volt;
using namespace Volt::Plugin;
using S = ReferenceCrystalMatchingService;

static const std::vector<OptionBinding<S>> bindings = {
    opt("--reference_crystal", "Perfect reference crystal (.lmp)", "", &S::setReferenceCrystal),
    opt("--anchor_species", "Species id whose dense sublattice anchors the grain frame", 1, &S::setAnchorSpecies),
    opt("--full_crystal_cutoff", "Reference bond cutoff in Angstrom (0 = auto-select)", 0.0, &S::setFullCrystalCutoff),
    opt("--topology_name", "Topology name written into the cluster table", "crystal", &S::setTopologyName),
    opt("--species_column", "Dump column carrying the atomic species", "species", &S::setSpeciesColumn),
    opt("--reference_topology", "Exotic YAML basename in --lattice_dir (anchor_species + reference_crystal)", "", &S::setReferenceTopology),
    opt("--lattice_dir", "Directory holding the exotic YAML + reference .lmp", "", &S::setLatticeDir),
};

VOLT_SERVICE_PLUGIN("volt-reference-crystal-matching", "Reference Crystal Matching", S, bindings)
