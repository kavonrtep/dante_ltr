# core_drapa test fixture

A 73 366 bp window of *Draparnaldia* (a green alga), used by
`tests/core.sh` to exercise core-domain mode on a genome **REXdb does
not cover at all**.

## Why this fixture exists

Lineage mode finds nothing here. On the full Drapa genome it reports
201 rank-`D` clusters and zero complete elements; on this window it
produces an empty GFF3:

```
$ ./dante_ltr -g dante.gff3 -s genome.fasta -o out
INFO: 100% of LTR protein domains do not reach lineage depth.
$ grep -vc '^#' out.gff3
0
```

Core mode seeds structurally, so the same input yields 4 ordered
Ty3/gypsy cores. That contrast is the point of the feature, and
`tests/core.sh` asserts it.

## Contents

```
genome.fasta   drapa_ctg1, 73 366 bp
dante.gff3     21 DANTE protein_domain features
```

Domain inventory: `RH 5, RT 4, INT 4, PROT 4, CHD 4` — four complete
gypsy cores, each with PROT 5' and CHD 3' of it, and **no GAG at all**.
That is the layout core mode is designed for: the accessory complement
is incomplete, so lineage mode's whole-element matching cannot fire,
while the ordered RT/RH/INT core is intact.

Every LTR domain here classifies no deeper than superfamily, which is
why `--fallback_mode` does not rescue it either — the domain complement
is missing, not merely under-resolved.

## Provenance

| | |
|---|---|
| assembly | `/mnt/ceph/454_data/Drapa/hifiasm/assembly_2025_07_30/DRA_2025_07_30/output/hifiasm_assembly.bp.p_ctg.gfa.fasta` |
| DANTE | `.../output/analysis/repeat_annotation/output/DANTE/DANTE.gff3` (raw, unfiltered) |
| region | `ptg000002l:19026193-19099558` |
| renamed to | `drapa_ctg1`, coordinates shifted to start at 1 |

The window was selected by scanning
`utils/measure_core_seeding.py --dump_seeds` output for the most compact
region containing ≥ 2 gypsy seeds that all carry a CHD within 2 kb 3' of
the core.

Regenerate with `./make_fixture.py` (sources are outside the repo, so
this is for maintenance, not CI).
