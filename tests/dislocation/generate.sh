set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"     
ROOT="$(cd "$REPO/.." && pwd)"                           
RCM="$REPO/build/Release/reference-crystal-matching"
DXA="$ROOT/opendxa/build/Release/opendxa"
LAT="$ROOT/opendxa/lattices"
FIX="$REPO/tests/fixtures"
OUT="$REPO/tests/dislocation"

# run_pipeline <name> <topo> <disloc_file> <ref_file> <extra_rcm_flags> [dxa_rescale]
#   Runs RCM (producer) on a dislocated dump against ref_file, then OpenDXA
#   (consumer) on the emitted contract. Writes all artifacts to OUT/<name>/.
#   dxa_rescale (optional) overrides the metric_rescale fed to OpenDXA — used for
#   A7/Bi where the affine isotropization corrupts Burgers and 1,1,1 is correct.
run_pipeline(){
    local name="$1" topo="$2" disloc="$3" ref="$4" extra="$5" dxa_rescale="${6:-}"
    local dir="$OUT/$name"
    local rescale
    rescale=$("$RCM" "$disloc" "$dir/${name}_rcm" \
        --reference_crystal "$ref" --anchor_species 1 \
        --species_column type --topology_name "$topo" $extra 2>&1 \
        | tee "$dir/${name}_rcm.log" \
        | sed -n 's/.*--metric_rescale \([0-9.,]*\).*/\1/p' | tail -1)
    local use_rescale="${dxa_rescale:-${rescale:-1,1,1}}"
    echo "  RCM rescale=${rescale:-1,1,1}  DXA rescale=$use_rescale"
    "$DXA" "$dir/${name}_rcm_annotated.dump" "$dir/${name}_dxa" \
        --clusters_table "$dir/${name}_rcm_clusters.table" \
        --clusters_transitions "$dir/${name}_rcm_cluster_transitions.table" \
        --neighbor_lattice "$dir/${name}_rcm_neighbor_lattice.parquet" \
        --reference_topology "$topo" --lattice_dir "$LAT" \
        --metric_rescale "$use_rescale" 2>&1 | tee "$dir/${name}_dxa.log" \
        | grep -iE "Found .* segments" || true
}

# ---------------------------------------------------------------------------
# Group 1 — multi-species crystals. ref == atomsk supercell source; edge_rm
# removes a half-plane (robust: no NaN). 5th field = extra RCM flags.
#   name : topology : lattice `a` (Burgers/cut) : supercell Nx Ny Nz : [flags]
MULTI_CASES=(
    "mgo:rocksalt:4.21:20 20 8:"
    "caf2:fluorite:5.46:16 16 6:"
    "srtio3:perovskite:3.905:16 16 6:"
    "cementite:cementite:5.089:14 10 12:"
    "cr3si:A15:4.56:16 16 8:"
    # C15 Laves: auto-cutoff isolates 2/3 of atoms; force 4.2 A to connect the
    # graph. NOTE: still 0 segments — RCM does not resolve the 24-site Laves
    # basis (see README). Kept as an honest negative case.
    "mgcu2:C15:7.05:8 8 10:--full_crystal_cutoff 4.2"
)

for entry in "${MULTI_CASES[@]}"; do
    IFS=':' read -r name topo a dup extra <<< "$entry"
    echo "========== $name ($topo) =========="
    dir="$OUT/$name"; mkdir -p "$dir"
    atomsk "$FIX/${name}_ref.lmp" -duplicate $dup "$dir/${name}_big.lmp" -ow >/dev/null 2>&1
    xhi=$(awk '/xhi/{print $2; exit}' "$dir/${name}_big.lmp")
    cx=$(python3 -c "print($xhi/2 + 0.03)")
    atomsk "$dir/${name}_big.lmp" -disloc "$cx" "$cx" edge_rm z y "$a" 0.3 \
        "$dir/${name}_disloc.lmp" -ow >/dev/null 2>&1
    rm -f "$dir/${name}_big.lmp"
    run_pipeline "$name" "$topo" "$dir/${name}_disloc.lmp" "$FIX/${name}_ref.lmp" "$extra"
done

# ---------------------------------------------------------------------------
# Group 2 — single-species Bravais lattices. The RCM reference is a *primitive*
# cell (committed in fixtures/, non-orthogonal for FCC/BCC/HCP/A7) so site
# assignment is unambiguous. The dislocation supercell is built from the
# conventional cell (atomsk --create) — or, for A7 which atomsk cannot create,
# from the committed hexagonal cell. HCP/A7 use screw dislocations because
# edge_rm crashes atomsk on hexagonal/rhombohedral boxes (NaN / End-of-file).
#
# A7 (arsenic/bismuth) is rhombohedral with a 2-atom primitive cell. Using the
# 6-atom hexagonal cell as the RCM reference collapses Burgers vectors to a full
# lattice vector (its sites are equivalent — same failure as a conventional FCC
# cell), so a7_prim.lmp (2-atom rhombohedral) is the reference, and DXA runs at
# metric_rescale 1,1,1 (the affine isotropization corrupts the Burgers here).
#   name : topology : source : dup : disloc burgers : ref_file : [flags] : [dxa_rescale]
#   source = "create:<atomsk create args>"  or  "ref:<fixture .lmp>"
BRAVAIS_CASES=(
    "fcc:fcc:create:fcc 3.615 Cu:20 20 8:edge_rm 3.615:fcc_prim.lmp::"
    "bcc:bcc:create:bcc 2.87 Fe:18 18 10:edge_rm 2.87:bcc_prim.lmp:--full_crystal_cutoff 3.0:"
    "hcp:hcp:create:hcp 3.21 5.21 Mg:16 16 6:screw 5.21:hcp_prim.lmp:--full_crystal_cutoff 4.0:"
    "a7:A7:ref:a7_ref.lmp:10 10 4:screw 10.55:a7_prim.lmp:--full_crystal_cutoff 4.0:1,1,1"
)

for entry in "${BRAVAIS_CASES[@]}"; do
    IFS=':' read -r name topo srckind srcspec dup discspec ref extra dxarescale <<< "$entry"
    echo "========== $name ($topo) =========="
    dir="$OUT/$name"; mkdir -p "$dir"
    # Build the conventional source cell.
    if [ "$srckind" = "create" ]; then
        atomsk --create $srcspec "$dir/${name}_conv.lmp" -ow >/dev/null 2>&1
        src="$dir/${name}_conv.lmp"
    else
        src="$FIX/$srcspec"
    fi
    atomsk "$src" -duplicate $dup "$dir/${name}_big.lmp" -ow >/dev/null 2>&1
    rm -f "$dir/${name}_conv.lmp"
    # Insert the dislocation. edge_rm closes cleanly at a centred core; screw is
    # placed off-centre (edge_rm crashes atomsk on hex/rhombohedral boxes).
    read -r dmode burgers <<< "$discspec"
    xhi=$(awk '/xhi/{print $2; exit}' "$dir/${name}_big.lmp")
    if [ "$dmode" = "screw" ]; then
        cx=$(python3 -c "print($xhi*0.37)"); cy=$(python3 -c "print($xhi*0.41)")
        atomsk "$dir/${name}_big.lmp" -disloc "$cx" "$cy" screw z y "$burgers" \
            "$dir/${name}_disloc.lmp" -ow >/dev/null 2>&1
    else
        cc=$(python3 -c "print($xhi/2 + 0.03)")
        atomsk "$dir/${name}_big.lmp" -disloc "$cc" "$cc" edge_rm z y "$burgers" 0.3 \
            "$dir/${name}_disloc.lmp" -ow >/dev/null 2>&1
    fi
    rm -f "$dir/${name}_big.lmp"
    run_pipeline "$name" "$topo" "$dir/${name}_disloc.lmp" "$FIX/$ref" "$extra" "$dxarescale"
done

echo "All artifacts under: $OUT"
