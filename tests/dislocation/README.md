# Dislocation end-to-end fixtures

Proof that the modularized **Reference Crystal Matching** producer + **OpenDXA**
consumer detect dislocations across a wide range of crystals — multi-species
oxides/intermetallics AND single-species Bravais lattices, including
non-orthogonal (hexagonal / rhombohedral) cells.

For each crystal, `generate.sh` inserts a dislocation with atomsk, runs it
through RCM → contract → OpenDXA, and stores every artifact here. Outputs are
committed so a run can be inspected without re-running atomsk/RCM/DXA.

## Two groups

**Multi-species** (group 1): the RCM reference is the same conventional cell
atomsk duplicates; an edge dislocation is inserted with `edge_rm` (removes a
half-plane — robust, no NaN) at a centred core.

**Single-species Bravais** (group 2): the RCM reference is a *primitive* cell
(committed `*_prim.lmp` / `a7_ref.lmp`, non-orthogonal for FCC/BCC/HCP/A7) so
basis-site assignment is unambiguous — a conventional FCC cell has 4 equivalent
sites that collapse the ideal vectors. The dislocation supercell is built from
the conventional cell (`atomsk --create`), except A7 which atomsk cannot create
(built from the committed reference). HCP/A7 use *screw* dislocations because
`edge_rm` crashes atomsk on hexagonal/rhombohedral boxes.

## Per-crystal artifacts (`<name>/`)

`<name>_disloc.lmp` (atomsk input), `<name>_rcm.log` + `<name>_rcm_*`
(producer contract: annotated dump, clusters/transitions tables,
neighbor_lattice, atoms.parquet), `<name>_dxa.log` + `<name>_dxa_*`
(consumer: dislocations, summary, defect_mesh parquet).

## Results

| Crystal | Structure | Symmetry | Species | metric_rescale | Segments | \|b\|dom |
|---------|-----------|----------|---------|----------------|----------|---------|
| mgo       | rocksalt   | cubic         | 2 | 1,1,1            | ~46 | 8.42 (=2a) |
| caf2      | fluorite   | cubic         | 2 | 1,1,1            | ~15 | 5.46 (=a) |
| srtio3    | perovskite | cubic         | 3 | 1,1,1            | 1   | 3.90 (=a) |
| cementite | cementite  | orthorhombic  | 2 | 0.755,1,0.671    | 1   | 6.74 (=b) |
| cr3si     | A15        | cubic         | 2 | 1,1,1            | 1   | 4.56 (=a) |
| fcc       | FCC        | cubic (prim.) | 1 | 1,0.866,0.816    | 1-2 | 3.63 (=a) |
| bcc       | BCC        | cubic (prim.) | 1 | 1,0.943,0.816    | 1   | 2.87 (=a) |
| hcp       | HCP        | hexagonal     | 1 | 0.616,0.534,1    | 1-2 | 5.21 (=c) |
| a7        | A7 arsenic | rhombohedral  | 1 | 1,1,1            | 1   | 10.55 (=c) |
| bi_np     | A7 bismuth (real MD NP, 221k atoms) | rhombohedral | 1 | 1,1,1 | ~221 | 4.55 (=a) |
| mgcu2     | C15 Laves  | cubic         | 2 | 1,1,1            | **0** ✗ | — |

Every working case recovers a **physical** Burgers vector equal to a lattice
parameter (the `|b|dom` column). A single straight line (1 segment) is the ideal
result for a clean dislocation; higher counts (mgo, bi_np) come from core
fragmentation / dissociation or, for the real NP, a genuine dislocation network —
with the dominant Burgers still a real lattice vector.

### A7 / rhombohedral — primitive cell is mandatory

A7 (arsenic, bismuth) MUST use a **2-atom primitive rhombohedral** reference
(`a7_prim.lmp`, `bi_prim.lmp`), not the 6-atom hexagonal cell. The hexagonal
cell has crystallographically equivalent sites that `assignBasisSites` cannot
distinguish — exactly the same failure as a conventional (4-site) FCC cell — so
the ideal vectors collapse to a full lattice vector and DXA returns a spurious
Burgers (|b|≈18 Å, ~4× the lattice parameter) on thousands of bogus segments.
With the primitive cell the Burgers come out physical (a7: |b|=c=10.55;
bi_np: |b|=a=4.55). A7 also runs at metric_rescale `1,1,1` — the affine
isotropization corrupts the Burgers for this low-symmetry cell.

The real Bi nanoparticle (`bi_np`, timestep 1.275M, 221021 atoms, free surfaces)
went from 9325 segments with garbage |b|=18.5 to **221 segments** with |b|=4.55
(=a) once the primitive reference replaced the hexagonal one.

### Known limitation — C15 Laves (mgcu2)

RCM does not resolve the C15 Laves phase: even in the defect-free bulk the
3-edge loop residual stays at a full lattice vector instead of ~0, so no Burgers
circuit closes (0 segments). Root cause is basis-site assignment for the 24-site
basis (Mg 8a + Cu 16d) — the per-atom shell match picks the wrong equivalent
site. Kept as an honest negative case. (Same class of bug as the A7 hexagonal
cell, but the C15 basis is not reducible to a single-orbit primitive cell, so
the primitive-cell fix does not apply.)

Segment counts vary slightly run-to-run (OpenDXA's TBB-parallel mesh tracing is
not bit-deterministic); the Burgers vectors do not.

## Regenerate

```bash
bash generate.sh   # needs atomsk on PATH + RCM/OpenDXA built (Release)
```
