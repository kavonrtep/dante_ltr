# element_at_sequence_end — regression fixture for issue #17

Unplaced contig `CAZZRL010000126.1` (47,569 bp) of *Matthiola sinuata*,
GCA_964656325.1 (ddMatSinu1.hap1.1), unchanged. `dante.gff3` holds the 5
DANTE domains of one Ty1/copia Angela element on `+` (GAG, PROT, INT, RT,
RH), in contig coordinates.

The RH domain ends on the last base of the contig, so `check_ranges()`
drops the only complete element. The empty list then reached
`get_ranges()`, which aborted with "non-numeric argument to binary
operator". The run must complete and report the element as rank D.
