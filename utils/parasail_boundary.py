"""Parasail-based anchored-extension LTR boundary refinement.

Library module — importable from refine_boundaries.py and any later tool.
Ports CARP's `global_local_aln.py` helpers (per-column scoring, optimal-
alignment extraction) and adds:

  * sequence-extraction geometry that respects + / - strand
  * TG/CA snap within a small window on the member's own raw sequence
  * a single `validate_motif()` choke point so switching policy
    (TG/CA -> none -> future per-lineage motif) is a one-line change

The module is a pure library: no I/O, no logging, no argparse.  The
caller is responsible for parsing GFF3, opening the genome, clustering,
and writing outputs.
"""

from __future__ import annotations

import re
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import parasail


# ============================================================
# Public data types
# ============================================================

@dataclass
class Member:
    """One LTR feature row from the input GFF3.

    Coordinates are 1-based inclusive (GFF3 convention).
    """
    feat_id: int
    chrom: str
    start: int
    end: int
    strand: str          # '+' / '-'
    role: str            # '5LTR' or '3LTR'
    parent: str          # parent transposable_element ID
    lineage: str         # Final_Classification


# ============================================================
# Sequence helpers
# ============================================================

_COMPLEMENT_TBL = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s: str) -> str:
    return s.translate(_COMPLEMENT_TBL)[::-1]


# ============================================================
# Motif validation — the single policy choke point
# ============================================================

# Supported motif policies in v1.  A future per-lineage column in
# lineage_domain_order.csv would resolve here as well; today the caller
# picks one policy globally.
SUPPORTED_MOTIFS = ("TG/CA", "none")


def validate_motif(boundary_pair: Optional[str], side: str, motif: str) -> Optional[bool]:
    """Validate a 2-bp boundary motif against a policy.

    boundary_pair : the 2-bp motif at the candidate boundary in
                    biological orientation (caller must already have
                    handled strand). May be None when out of range.
    side          : "5" or "3"
    motif         : policy. One of SUPPORTED_MOTIFS.

    Returns:
        True / False                 — TG/CA policy, motif present / absent
        True                         — motif=none (always validates)
        None                         — boundary_pair is None / cannot evaluate
    """
    if motif == "none":
        return True
    if boundary_pair is None:
        return None
    if motif == "TG/CA":
        target = "TG" if side == "5" else "CA"
        return boundary_pair == target
    raise ValueError(f"Unsupported boundary motif: {motif!r}")


# ============================================================
# Parasail anchored-extension primitives (ported from CARP)
# ============================================================

def _pick_func(end_anchor: str):
    """Return parasail trace function for the given end-fixed setup.

    end_anchor "5": fix 5' on both, free 3' (sg_qe_de_*)
    end_anchor "3": fix 3' on both, free 5' (sg_qb_db_*)
    """
    if end_anchor == "5":
        return parasail.sg_qe_de_trace_striped_16
    return parasail.sg_qb_db_trace_striped_16


def _make_matrix(match: int, mismatch: int):
    return parasail.matrix_create("ACGT", match, mismatch)


def _per_column_scores(result, match=2, mismatch=-2,
                       gap_open=12, gap_extend=3):
    """Affine-gap per-column scoring on a parasail traceback comp string."""
    tb = result.traceback
    comp = tb.comp
    q, r = tb.query, tb.ref
    n = len(comp)
    scores = [0] * n
    in_gap_q = False
    in_gap_r = False
    for i in range(n):
        c = comp[i]; qc = q[i]; rc = r[i]
        if c == '|':
            in_gap_q = False; in_gap_r = False
            scores[i] = match
        elif c == '.':
            in_gap_q = False; in_gap_r = False
            scores[i] = mismatch
        elif c == ' ':
            if qc == '-' and rc != '-':
                if not in_gap_q:
                    scores[i] = -(gap_open + gap_extend); in_gap_q = True
                else:
                    scores[i] = -gap_extend
                in_gap_r = False
            elif qc != '-' and rc == '-':
                if not in_gap_r:
                    scores[i] = -(gap_open + gap_extend); in_gap_r = True
                else:
                    scores[i] = -gap_extend
                in_gap_q = False
    return scores


def _extract_optimal(result, scores, end="5"):
    """Walk cumulative score from the fixed end; report degapped lengths to peak."""
    tb = result.traceback
    q, r = tb.query, tb.ref
    n = len(scores)
    if end == "5":
        cs = 0; best = -10**9; bp = 0
        for i, s in enumerate(scores):
            cs += s
            if cs > best:
                best = cs; bp = i
        trimmed_q = q[:bp + 1]; trimmed_r = r[:bp + 1]
    else:
        cs = 0; best = -10**9; bp = n - 1
        for i in range(n - 1, -1, -1):
            cs += scores[i]
            if cs > best:
                best = cs; bp = i
        trimmed_q = q[bp:]; trimmed_r = r[bp:]
    dq = trimmed_q.replace('-', '')
    dr = trimmed_r.replace('-', '')
    return {"max_score": best, "max_pos": bp,
            "degapped_query_len": len(dq), "degapped_ref_len": len(dr)}


def align_pair(seq1: str, seq2: str, end: str,
               match: int = 2, mismatch: int = -2,
               gap_open: int = 12, gap_extend: int = 3) -> Dict:
    """Run parasail semi-global alignment with one end fixed; return summary.

    Returns dict with keys: max_score, max_pos, degapped_query_len,
    degapped_ref_len, parasail_score.
    """
    func = _pick_func(end)
    mat = _make_matrix(match, mismatch)
    res = func(seq1, seq2, gap_open, gap_extend, mat)
    scores = _per_column_scores(res, match=match, mismatch=mismatch,
                                gap_open=gap_open, gap_extend=gap_extend)
    out = _extract_optimal(res, scores, end=end)
    out["parasail_score"] = res.score
    return out


# ============================================================
# Boundary geometry on a member's raw genomic sequence
# ============================================================

def get_boundary_pair(genome, chrom: str, strand: str, pos: int,
                      side: str, chrom_len: int) -> Optional[str]:
    """Return the 2-bp boundary motif in biological orientation at the
    candidate corrected boundary `pos` (1-based genomic), or None if
    out of range.

      side="5", + strand : motif at genome[pos .. pos+1]
      side="5", - strand : revcomp(genome[pos-1 .. pos])
      side="3", + strand : motif at genome[pos-1 .. pos]
      side="3", - strand : revcomp(genome[pos .. pos+1])
    """
    try:
        if side == "5":
            if strand == '+':
                if pos < 1 or pos + 1 > chrom_len:
                    return None
                return str(genome[chrom][pos - 1:pos + 1]).upper()
            else:
                if pos - 1 < 1 or pos > chrom_len:
                    return None
                return revcomp(str(genome[chrom][pos - 2:pos]).upper())
        else:  # side == "3"
            if strand == '+':
                if pos - 1 < 1 or pos > chrom_len:
                    return None
                return str(genome[chrom][pos - 2:pos]).upper()
            else:
                if pos < 1 or pos + 1 > chrom_len:
                    return None
                return revcomp(str(genome[chrom][pos - 1:pos + 1]).upper())
    except Exception:
        return None


def snap_to_motif(genome, m: Member, corrected_g: Optional[int], side: str,
                  chrom_len: int, window: int, motif: str
                  ) -> Tuple[Optional[int], int, Optional[bool]]:
    """Search corrected_g and ± `window` bp for nearest position whose
    boundary motif matches the active policy.  Closest-to-original wins
    (offset 0 first, then expanding alternately to ±1, ±2, ...).

    Returns (snapped_pos, offset, motif_match) where:
      snapped_pos : the position whose boundary motif matches policy,
                    or `corrected_g` if no match in window.
      offset      : snapped_pos - corrected_g
      motif_match : True / False / None (None when motif=none and there
                    is no policy to validate against -- callers should
                    treat that as "validated by parasail evidence alone")

    If motif=="none", we return (corrected_g, 0, True) without searching.
    """
    if corrected_g is None:
        return None, 0, False
    if motif == "none":
        return corrected_g, 0, True

    target = "TG" if side == "5" else "CA"
    pair0 = get_boundary_pair(genome, m.chrom, m.strand, corrected_g,
                              side, chrom_len)
    if pair0 == target:
        return corrected_g, 0, True
    if window <= 0:
        return corrected_g, 0, False
    for d in range(1, window + 1):
        for sgn in (-1, 1):
            cand = corrected_g + sgn * d
            pair = get_boundary_pair(genome, m.chrom, m.strand, cand,
                                     side, chrom_len)
            if pair == target:
                return cand, sgn * d, True
    return corrected_g, 0, False


# ============================================================
# Per-member sequence extraction for parasail
# ============================================================

def extract_seq_for_side(member: Member, genome, side: str,
                         anchor_len: int, flank_len: int,
                         genome_lens: Dict[str, int]) -> Dict:
    """Return per-member input for parasail anchored extension on `side`.

    Layout in biological orientation:

      side="5":  [flank_len bp upstream genomic] [anchor at LTR 5' end]
                 (anchor is at the 3' end of the returned string -> end="3")

      side="3":  [anchor at LTR 3' end] [flank_len bp downstream genomic]
                 (anchor is at the 5' end of the returned string -> end="5")

    Caller should call align_pair(seq_i, seq_j, end="3" if side=="5" else "5").
    """
    glen = genome_lens[member.chrom]
    if member.strand == '+':
        if side == "5":
            anchor_g_start = member.start
            anchor_g_end = min(member.end, member.start + anchor_len - 1)
            flank_g_start = max(1, member.start - flank_len)
            flank_g_end = member.start - 1
            anchor_seq = str(genome[member.chrom][anchor_g_start - 1:anchor_g_end]).upper()
            flank_seq = ("" if flank_g_end < flank_g_start else
                         str(genome[member.chrom][flank_g_start - 1:flank_g_end]).upper())
            return {"seq": flank_seq + anchor_seq,
                    "anchor": len(anchor_seq), "flank": len(flank_seq),
                    "boundary_g": member.start}
        else:
            anchor_g_start = max(member.start, member.end - anchor_len + 1)
            anchor_g_end = member.end
            flank_g_start = member.end + 1
            flank_g_end = min(glen, member.end + flank_len)
            anchor_seq = str(genome[member.chrom][anchor_g_start - 1:anchor_g_end]).upper()
            flank_seq = ("" if flank_g_end < flank_g_start else
                         str(genome[member.chrom][flank_g_start - 1:flank_g_end]).upper())
            return {"seq": anchor_seq + flank_seq,
                    "anchor": len(anchor_seq), "flank": len(flank_seq),
                    "boundary_g": member.end}
    else:  # '-' strand
        if side == "5":
            anchor_g_start = max(member.start, member.end - anchor_len + 1)
            anchor_g_end = member.end
            flank_g_start = member.end + 1
            flank_g_end = min(glen, member.end + flank_len)
            anchor_seq = str(genome[member.chrom][anchor_g_start - 1:anchor_g_end]).upper()
            flank_seq = ("" if flank_g_end < flank_g_start else
                         str(genome[member.chrom][flank_g_start - 1:flank_g_end]).upper())
            return {"seq": revcomp(anchor_seq + flank_seq),
                    "anchor": len(anchor_seq), "flank": len(flank_seq),
                    "boundary_g": member.end}
        else:
            anchor_g_start = member.start
            anchor_g_end = min(member.end, member.start + anchor_len - 1)
            flank_g_start = max(1, member.start - flank_len)
            flank_g_end = member.start - 1
            anchor_seq = str(genome[member.chrom][anchor_g_start - 1:anchor_g_end]).upper()
            flank_seq = ("" if flank_g_end < flank_g_start else
                         str(genome[member.chrom][flank_g_start - 1:flank_g_end]).upper())
            return {"seq": revcomp(flank_seq + anchor_seq),
                    "anchor": len(anchor_seq), "flank": len(flank_seq),
                    "boundary_g": member.start}


# ============================================================
# Per-cluster pairwise parasail
# ============================================================

def extract_pool_for_side(members: List[Member], role_keep: str, side: str,
                          genome, genome_lens: Dict[str, int],
                          anchor_len: int, flank_len: int) -> List[Dict]:
    """Build a list of per-member info dicts (seq + anchor metadata) for
    members whose `.role` matches `role_keep`.

    The geometry returned by extract_seq_for_side() depends on the
    member's coordinates and strand only — not the role label — so this
    function can be used to construct *opposite-role* anchor pools (the
    inner-edge pool from §2 of the v2 design): pass `role_keep` as the
    role you want to keep (e.g. '3LTR' when refining 5'LTR boundaries
    via the inner pool) and the same `side` value the query is using.
    """
    out: List[Dict] = []
    for m in members:
        if m.role != role_keep:
            continue
        info = extract_seq_for_side(m, genome, side, anchor_len, flank_len,
                                    genome_lens)
        info["member"] = m
        out.append(info)
    return out


def project_corrected_g(member: Member, side: str,
                        ext_len: Optional[int],
                        anchor_used: int) -> Optional[int]:
    """Map a per-query degapped extension length to a 1-based genomic
    coordinate for the refined boundary.  Returns None if ext_len is None.
    """
    if ext_len is None:
        return None
    flank_ext = ext_len - anchor_used
    if side == "5":
        return (member.start - flank_ext) if member.strand == "+" \
               else (member.end + flank_ext)
    return (member.end + flank_ext) if member.strand == "+" \
           else (member.start - flank_ext)


def aggregate_extension(query_info: Dict, anchor_infos: List[Dict],
                        side: str,
                        match: int = 2, mismatch: int = -2,
                        gap_open: int = 12, gap_extend: int = 3,
                        score_threshold: int = 20,
                        min_num_alignments: int = 3,
                        max_pairs: int = 0) -> Dict:
    """Pairwise-align a single query against each anchor; aggregate the
    Nth-largest degapped query extension as the consensus call.

    Returns a dict with keys:
      ext_len     : Nth-largest extension (or None if n_pairs < N)
      n_pairs     : number of alignments that passed score_threshold
      score_max   : max parasail score across pairs (None if no pairs)
    """
    end_param = "3" if side == "5" else "5"
    if (query_info["anchor"] < 10 or
            len(query_info["seq"]) < query_info["anchor"] + 5):
        return {"ext_len": None, "n_pairs": 0, "score_max": None}

    # Subsample if too many anchors
    anchors = anchor_infos
    if max_pairs and len(anchors) > max_pairs:
        import random
        rng = random.Random(123)
        anchors = rng.sample(anchors, max_pairs)

    ext_lens: List[int] = []
    score_max: Optional[int] = None
    for a in anchors:
        if a["anchor"] < 10 or len(a["seq"]) < a["anchor"] + 5:
            continue
        # Skip self-alignment if query and anchor are the same member
        if a.get("member") is not None and \
                a["member"].feat_id == query_info.get("member") and \
                query_info.get("member") is not None and \
                a["member"].feat_id == query_info["member"].feat_id:
            continue
        try:
            r = align_pair(query_info["seq"], a["seq"], end_param,
                           match=match, mismatch=mismatch,
                           gap_open=gap_open, gap_extend=gap_extend)
        except Exception:
            continue
        if r["max_score"] < score_threshold:
            continue
        ext_lens.append(r["degapped_query_len"])
        if score_max is None or r["max_score"] > score_max:
            score_max = r["max_score"]

    if len(ext_lens) < min_num_alignments:
        return {"ext_len": None, "n_pairs": len(ext_lens), "score_max": score_max}
    sel = sorted(ext_lens, reverse=True)[min_num_alignments - 1]
    return {"ext_len": sel, "n_pairs": len(ext_lens), "score_max": score_max}


def process_cluster_side(cluster_members: List[Member], side: str,
                         genome, genome_lens: Dict[str, int],
                         anchor_len: int, flank_len: int,
                         match: int = 2, mismatch: int = -2,
                         gap_open: int = 12, gap_extend: int = 3,
                         score_threshold: int = 20,
                         min_num_alignments: int = 3,
                         snap_window: int = 5,
                         motif: str = "TG/CA",
                         max_pairs: int = 0
                         ) -> Tuple[List[Dict], int]:
    """Run pairwise parasail across the cluster members of the role
    matching `side` (5'LTRs for side='5', 3'LTRs for side='3') and
    aggregate per-member extension lengths.

    Returns (per_member_rows, n_pairs_used).  Each row is a dict; see
    keys at the bottom of this function.

    Mixing roles is incorrect *for v1*: the same-role pool requires both
    sides of the alignment have variable genomic flanks.  v2 uses a
    different geometry for opposite-role pools — see
    `extract_pool_for_side` and `aggregate_extension` for the reusable
    primitives.  Callers must pass cluster_members containing both 5'
    and 3'LTRs; this function filters to the appropriate role per side.
    """
    role_keep = "5LTR" if side == "5" else "3LTR"
    seqinfo = extract_pool_for_side(cluster_members, role_keep, side,
                                    genome, genome_lens, anchor_len, flank_len)

    # Pairwise all-vs-all uses the same engine via aggregate_extension,
    # but we need each pair to contribute to BOTH endpoints, so we keep
    # the explicit pair loop here for the v1 path.
    n = len(seqinfo)
    pair_indices = [(i, j) for i in range(n) for j in range(i + 1, n)
                    if seqinfo[i]["anchor"] >= 10 and seqinfo[j]["anchor"] >= 10
                    and len(seqinfo[i]["seq"]) >= anchor_len + 5
                    and len(seqinfo[j]["seq"]) >= anchor_len + 5]
    if max_pairs and len(pair_indices) > max_pairs:
        import random
        rng = random.Random(123)
        pair_indices = rng.sample(pair_indices, max_pairs)

    end_param = "3" if side == "5" else "5"
    pool: Dict[int, List[int]] = defaultdict(list)
    n_pairs_used = 0
    for i, j in pair_indices:
        s1 = seqinfo[i]["seq"]; s2 = seqinfo[j]["seq"]
        try:
            r = align_pair(s1, s2, end_param,
                           match=match, mismatch=mismatch,
                           gap_open=gap_open, gap_extend=gap_extend)
        except Exception:
            continue
        if r["max_score"] < score_threshold:
            continue
        pool[i].append(r["degapped_query_len"])
        pool[j].append(r["degapped_ref_len"])
        n_pairs_used += 1

    out: List[Dict] = []
    for k in range(n):
        m = seqinfo[k]["member"]; anchor_used = seqinfo[k]["anchor"]
        lengths = pool.get(k, [])
        if len(lengths) >= min_num_alignments:
            sel = sorted(lengths, reverse=True)[min_num_alignments - 1]
            ext_len = sel
        else:
            sel = None
            ext_len = None
        corrected_g = project_corrected_g(m, side, ext_len, anchor_used)

        chrom_len = genome_lens[m.chrom]
        motif_at_parasail = None
        motif_ok = None
        snap_offset = 0
        corrected_g_final = corrected_g
        if corrected_g is not None:
            pair_at_parasail = get_boundary_pair(
                genome, m.chrom, m.strand, corrected_g, side, chrom_len)
            motif_at_parasail = validate_motif(pair_at_parasail, side, motif)
            corrected_g_final, snap_offset, motif_ok = snap_to_motif(
                genome, m, corrected_g, side, chrom_len, snap_window, motif)

        out.append({
            "member_idx_in_cluster": k,
            "member_feat_id": m.feat_id,
            "chrom": m.chrom, "start": m.start, "end": m.end,
            "strand": m.strand, "role": m.role,
            "side": side,
            "n_pairs": len(lengths),
            "anchor_used": anchor_used,
            "ext_len_nth": sel,
            "flank_extension": (ext_len - anchor_used) if ext_len is not None else None,
            "annotated_g": (m.start if (side == "5" and m.strand == '+') or
                                       (side == "3" and m.strand == '-')
                            else m.end),
            "corrected_g_parasail": corrected_g,
            "motif_at_parasail": motif_at_parasail,
            "snap_offset": snap_offset,
            "corrected_g": corrected_g_final,
            "motif_ok": motif_ok,
        })
    return out, n_pairs_used
