#!/usr/bin/env bash
# Generate dislocation end-to-end fixtures: insert an edge dislocation into each
# reference crystal with atomsk, run it through RCM (the producer) and then
# OpenDXA (the consumer), and keep every artifact under tests/dislocation/<name>/.
#
# This proves the modularized RCM + OpenDXA detect dislocations in crystals that
# are NOT forsterite. Outputs are committed so the run can be inspected without
# re-running atomsk/RCM/DXA. Re-run this script to regenerate.
#
# Requires: atomsk on PATH, RCM + OpenDXA already built (Release).
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"          # reference-crystal-matching/
ROOT="$(cd "$REPO/.." && pwd)"                                       # repos-git/
RCM="$REPO/build/Release/reference-crystal-matching"
DXA="$ROOT/opendxa/build/Release/opendxa"
LAT="$ROOT/opendxa/lattices"
FIX="$REPO/tests/fixtures"
OUT="$REPO/tests/dislocation"

# name : topology : lattice param `a` (Burgers length / cut) : supercell Nx Ny Nz
#   line along z, glide cut along y, edge_rm removes a half-plane (robust: no NaN).
# name : topology : lattice param `a` (Burgers length / cut) : supercell Nx Ny Nz : [extra RCM flags]
#   line along z, glide cut along y, edge_rm removes a half-plane (robust: no NaN).
#   The optional 5th field passes extra flags to RCM (e.g. a forced cutoff for
#   dense-basis lattices whose auto-selected cutoff isolates atoms).
CASES=(
    "mgo:rocksalt:4.21:20 20 8:"
    "caf2:fluorite:5.46:16 16 6:"
    "srtio3:perovskite:3.905:16 16 6:"
    "cementite:cementite:5.089:14 10 12:"
    "cr3si:A15:4.56:16 16 8:"
    # C15 Laves: auto-cutoff isolates 2/3 of atoms; force 4.2 A to connect the
    # graph. NOTE: still yields 0 segments — RCM does not resolve the 24-site
    # Laves basis (see README). Kept as an honest negative case.
    "mgcu2:C15:7.05:8 8 10:--full_crystal_cutoff 4.2"
)

for entry in "${CASES[@]}"; do
    IFS=':' read -r name topo a dup extra <<< "$entry"
    echo "========== $name ($topo) =========="
    dir="$OUT/$name"
    mkdir -p "$dir"

    # 1) Big supercell + edge dislocation (half-plane removed). -ow overwrites
    #    prior outputs so re-runs don't block on atomsk's confirm prompt.
    atomsk "$FIX/${name}_ref.lmp" -duplicate $dup "$dir/${name}_big.lmp" -ow >/dev/null 2>&1
    xhi=$(awk '/xhi/{print $2; exit}' "$dir/${name}_big.lmp")
    cx=$(python3 -c "print($xhi/2 + 0.03)")
    atomsk "$dir/${name}_big.lmp" -disloc "$cx" "$cx" edge_rm z y "$a" 0.3 \
        "$dir/${name}_disloc.lmp" -ow >/dev/null 2>&1
    rm -f "$dir/${name}_big.lmp"

    # 2) RCM producer -> contract. Capture the reported metric_rescale.
    rescale=$("$RCM" "$dir/${name}_disloc.lmp" "$dir/${name}_rcm" \
        --reference_crystal "$FIX/${name}_ref.lmp" --anchor_species 1 \
        --species_column type --topology_name "$topo" $extra 2>&1 \
        | tee "$dir/${name}_rcm.log" \
        | sed -n 's/.*--metric_rescale \([0-9.,]*\).*/\1/p' | tail -1)
    echo "  RCM metric_rescale=$rescale"

    # 3) OpenDXA consumer.
    "$DXA" "$dir/${name}_rcm_annotated.dump" "$dir/${name}_dxa" \
        --clusters_table "$dir/${name}_rcm_clusters.table" \
        --clusters_transitions "$dir/${name}_rcm_cluster_transitions.table" \
        --neighbor_lattice "$dir/${name}_rcm_neighbor_lattice.parquet" \
        --reference_topology "$topo" --lattice_dir "$LAT" \
        --metric_rescale "$rescale" 2>&1 | tee "$dir/${name}_dxa.log" \
        | grep -iE "Found .* segments" || true
done

echo "All artifacts under: $OUT"
