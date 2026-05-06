"""TSD (Target Site Duplication) re-detection at refined LTR boundaries.

Pure-Python port of the algorithm in `utils/ltr_utils.R:651-707`
(`evaluate_ltr` TSD branch).  Validated bit-for-bit against R on both the
`at` (Alyr) and `g2` validation datasets — see
`docs/refine_v2_analysis.md` §4.4 (1 026 / 1 920 elements with TSD all
match exact-sequence; 0 mismatches; 4 / 1 884 fuzzy edge cases on g2).

Public entry point:

    detect_tsd(genome, chrom, ltr5_outer, ltr3_outer, strand,
               max_len=8, min_exact=4, min_fuzzy=5)
        -> (length, sequence, fuzzy_flag)

Returns (0, "", False) when no TSD is found within the policy.
"""
from __future__ import annotations

from typing import Tuple


_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def detect_tsd(genome, chrom: str, ltr5_outer: int, ltr3_outer: int,
               strand: str, max_len: int = 8,
               min_exact: int = 4, min_fuzzy: int = 5
               ) -> Tuple[int, str, bool]:
    """Detect TSD flanking the (refined) outer boundaries of the two LTRs.

    Parameters
    ----------
    genome : pyfaidx.Fasta or compatible (genome[chrom][a:b] returns str)
    chrom  : chromosome / contig name
    ltr5_outer : 1-based genomic coord of the 5'LTR's outer 5' bp (the TG)
    ltr3_outer : 1-based genomic coord of the 3'LTR's outer 3' bp (the CA)
    strand     : '+' or '-' (the element's strand)
    max_len    : largest TSD length considered (default 8)
    min_exact  : minimum length for an exact-match TSD (default 4)
    min_fuzzy  : minimum length for a 1-mismatch fuzzy match (default 5)

    Returns
    -------
    (length, sequence, fuzzy_flag)
        length = 0 when no TSD found.
        sequence is in biological orientation (TSD_L's bases).  When
        fuzzy_flag is True and the two halves differ, sequence is just
        the L-side bases (callers needing both halves should re-extract).
    """
    chrom_len = len(genome[chrom])
    if strand == "+":
        # bp immediately upstream of 5'LTR start (TSD_L) and downstream
        # of 3'LTR end (TSD_R) in genomic forward orientation.
        l_start, l_end = ltr5_outer - max_len, ltr5_outer - 1
        r_start, r_end = ltr3_outer + 1, ltr3_outer + max_len
        if l_start < 1 or l_end > chrom_len or r_start < 1 or r_end > chrom_len:
            return (0, "", False)
        l_seq = str(genome[chrom][l_start - 1:l_end]).upper()
        r_seq = str(genome[chrom][r_start - 1:r_end]).upper()
    else:
        # On - strand the element reads 5'->3' from genomic high to low.
        # TSD_L (biological upstream of element) sits at + strand right
        # of ltr5_outer; TSD_R (biological downstream) sits at + strand
        # left of ltr3_outer.  Reverse-complement so both halves are in
        # biological orientation with TSD_L's rightmost bp closest to
        # the LTR start.
        l_start, l_end = ltr5_outer + 1, ltr5_outer + max_len
        r_start, r_end = ltr3_outer - max_len, ltr3_outer - 1
        if l_start < 1 or l_end > chrom_len or r_start < 1 or r_end > chrom_len:
            return (0, "", False)
        l_seq = _revcomp(str(genome[chrom][l_start - 1:l_end]).upper())
        r_seq = _revcomp(str(genome[chrom][r_start - 1:r_end]).upper())

    # Longest exact match wins.  R also scans length 1..3 then masks
    # them; we just start from min_exact.
    for n in range(max_len, min_exact - 1, -1):
        if l_seq[-n:] == r_seq[:n]:
            return (n, l_seq[-n:], False)

    # Fallback: ≤1 mismatch, length ≥ min_fuzzy.
    for n in range(max_len, min_fuzzy - 1, -1):
        l = l_seq[-n:]; r = r_seq[:n]
        if sum(1 for a, b in zip(l, r) if a != b) <= 1:
            return (n, l, True)

    return (0, "", False)


def _selftest() -> None:
    """Tiny self-test: synthetic chromosome with a known TSD."""
    class _FakeChrom:
        def __init__(self, s): self._s = s
        def __getitem__(self, sl): return self._s[sl]
        def __len__(self): return len(self._s)

    class _FakeFasta(dict):
        pass

    # Construct: ...AAAAA-CGTACGT[TG...CA]ACTGACTG[TG...CA]CGTACGT-AAAAA...
    # TSD = "CGTACGT" (7 bp) flanking both LTRs on + strand.
    upstream = "A" * 100 + "CGTACGT"           # bp 1..107
    ltr5     = "TG" + "X" * 100 + "CA"         # bp 108..211 (104 bp LTR)
    interior = "ACTGACTG" * 50                  # bp 212..611
    ltr3     = "TG" + "Y" * 100 + "CA"         # bp 612..715
    downstream = "CGTACGT" + "A" * 100         # bp 716..822
    seq = (upstream + ltr5 + interior + ltr3 + downstream).replace("X", "T").replace("Y", "T")

    g = _FakeFasta()
    g["chr1"] = _FakeChrom(seq)

    # 5'LTR outer = position 108 (first 'T' of TG); 3'LTR outer = 715 (last 'A' of CA)
    n, s, fz = detect_tsd(g, "chr1", 108, 715, "+")
    assert (n, s, fz) == (7, "CGTACGT", False), (n, s, fz)

    # - strand mirror: build an element where the same TSD flanks but the
    # whole thing is on the - strand. ltr5_outer would then be the genomic
    # high bp of the 5'LTR.  Construct by reverse-complementing the seq.
    seq_rc = _revcomp(seq)
    g["chr2"] = _FakeChrom(seq_rc)
    rc_len = len(seq_rc)
    # On the rc, the element's 5'LTR outer (biological 5' end) = rc_len - 108 + 1
    # 3'LTR outer = rc_len - 715 + 1
    ltr5_outer_rc = rc_len - 108 + 1
    ltr3_outer_rc = rc_len - 715 + 1
    n, s, fz = detect_tsd(g, "chr2", ltr5_outer_rc, ltr3_outer_rc, "-")
    assert (n, s, fz) == (7, "CGTACGT", False), (n, s, fz)

    print("tsd_redetect self-test OK")


if __name__ == "__main__":
    _selftest()
