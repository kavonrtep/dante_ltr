# core_weak — regression fixture for issue #14

`dante.gff3` is `tests/data/smoke/dante.gff3` with the element's five
domains (RH, RT, INT, PROT, GAG) weakened to `Similarity=0.6`,
`Relat_Length=0.5`. They pass the seed filter
(`--min_relative_length_core` 0.3) but fail the block filter
(`--min_relative_length` 0.6). The first domain, at 8 kb, is kept strong
and copied 1500 bp downstream so the block set has two domains outside
the element.

Use with `tests/data/smoke/genome.fasta`. Before the fix, core mode
re-collected the element's domains from the block set only, found none,
and aborted in `get_te_gff3()`. The element must now be reported with
its RH/RT/INT core as `protein_domain` children.

`dante_single_block.gff3` is the same without the copied domain, so
exactly one domain passes the block filter. `get_domain_clusters_alt()`
used to fail on a one-domain input when building the rank-D track.
