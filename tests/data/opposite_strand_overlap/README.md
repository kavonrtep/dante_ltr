# opposite_strand_overlap — regression fixture for issue #15

30 kb of *Prunus brigantina* chromosome 7, GCA_964267205.1 (drPruBrig1.1),
`OZ184548.1:21955001-21985000`, renamed `contig1`. `dante.gff3` holds the
7 DANTE 0.2.11 domains in that window, shifted to window coordinates.

A Ty1/copia Ivana element on `+` (GAG, PROT, INT, RT; complete only with
`-M 1`) has a `-` strand Ikeros RH starting inside its RT. Neighbour
lookup is strand-blind, so the right-hand BLAST window got a negative
offset and `get_ranges_right()` aborted with "each range must have an end
that is greater or equal to its start minus one". The run must complete.
It also yields no element of rank DL or higher, which used to break
`dante_ltr_summary`.
