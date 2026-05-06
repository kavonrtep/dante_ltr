"""Per-element LTR boundary refinement engine (v2).

Reads a DANTE_LTR GFF3 + reference FASTA, runs MMseqs2 per-lineage
clustering, refines each cluster member's outer-facing LTR boundary
using the *inner-primary policy*, applies a per-element TSD-loss
revert gate, and emits:

  * <prefix>_refined.gff3            -- refined annotations + TSD child rows
  * <prefix>_per_element.tsv         -- per-element refinement record
  * <prefix>_clusters.tsv            -- cluster manifest with stats
  * <prefix>_run.json                -- parameters + timing

Design background and term definitions:
  docs/refine_v2_analysis.md
  docs/refine_v2_implementation_plan.md

Pipeline (per cluster):

  1. For each member × side (5'/3' outer boundary), run parasail
     anchored extension TWICE: once with the OUTER-EDGE pool
     (same-role anchors) and once with the INNER-EDGE pool
     (opposite-role anchors used at their inner edge).
  2. Inner-primary policy: take the inner-pool coordinate when its
     TG/CA motif validates; else keep DANTE_LTR's original.
  3. Confidence label per side: dual / divergent / inner_only / unrefined.
  4. Per-element TSD recheck at proposed coordinates.
  5. If applying the proposals would *destroy* the originally-detected
     TSD on this element, revert all changed sides on the element.
  6. Optional MAFFT change-point fallback for clusters whose inner-pool
     validation rate is low; MAFFT proposals are subject to the same
     motif + TSD gates.

The MAFFT fallback is implemented as an out-of-process Rscript call to
utils/refine_mafft_fallback.R; this Python module is the orchestrator.

This module is invoked through the dante_ltr_refine CLI wrapper.
"""

from __future__ import annotations

import argparse
import gzip
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

# Local lib
HERE = os.path.dirname(os.path.realpath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from parasail_boundary import (  # type: ignore
    Member,
    aggregate_extension,
    extract_pool_for_side,
    get_boundary_pair,
    project_corrected_g,
    revcomp,
    snap_to_motif,
    validate_motif,
)
from tsd_redetect import detect_tsd  # type: ignore


logger = logging.getLogger("dante_ltr_refine")


# ============================================================
# GFF3 parsing
# ============================================================

_ATTR_RE = re.compile(r'([A-Za-z0-9_]+)=([^;]+)')


def _attr(attr_str: str, key: str) -> Optional[str]:
    for k, v in _ATTR_RE.findall(attr_str):
        if k == key:
            return v
    return None


def _open_gff(path: str):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def parse_gff3_for_refinement(path: str
                              ) -> Tuple[List[Member], Dict[str, Dict],
                                         List[str]]:
    """Parse a DANTE_LTR GFF3.

    Returns:
      members  : list of Member (one per long_terminal_repeat row)
      te_meta  : dict keyed by TE id -> {classification, rank,
                                          tsd_seq, tsd_len}
                 tsd_seq is the input GFF's TSD= attribute (or "")
                 tsd_len is parsed length (0 when TSD=not_found)
      raw_lines: list of raw lines (full file) to be replayed when emitting
                 the refined GFF3.  Header lines and unrelated rows are
                 passed through verbatim.
    """
    members: List[Member] = []
    te_meta: Dict[str, Dict] = {}
    raw_lines: List[str] = []

    feat_idx = 0
    with _open_gff(path) as f:
        for line in f:
            raw_lines.append(line)
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                continue
            parts = stripped.split('\t')
            if len(parts) < 9:
                continue
            seqname, _src, ftype, s1, e1, _sc, strand, _ph, attrs = parts
            if ftype == 'transposable_element':
                te_id = _attr(attrs, 'ID')
                if te_id:
                    tsd_str = (_attr(attrs, 'TSD') or "").strip()
                    if tsd_str in ("", "not_found"):
                        tsd_seq = ""; tsd_len = 0
                    elif "/" in tsd_str:
                        # fuzzy: stored as "L/R" — len is bp count of one side
                        tsd_seq = tsd_str
                        tsd_len = len(tsd_str.split("/", 1)[0])
                    else:
                        tsd_seq = tsd_str
                        tsd_len = len(tsd_str)
                    te_meta[te_id] = {
                        'classification': _attr(attrs, 'Final_Classification') or "",
                        'rank': _attr(attrs, 'Rank') or "",
                        'tsd_seq': tsd_seq,
                        'tsd_len': tsd_len,
                    }
            elif ftype == 'long_terminal_repeat':
                role = _attr(attrs, 'LTR') or ''
                parent = _attr(attrs, 'Parent') or ''
                if role not in ('5LTR', '3LTR'):
                    continue
                members.append(Member(
                    feat_id=feat_idx,
                    chrom=seqname,
                    start=int(s1), end=int(e1),
                    strand=strand if strand in ('+', '-') else '+',
                    role=role, parent=parent, lineage="",
                ))
                feat_idx += 1

    for m in members:
        m.lineage = te_meta.get(m.parent, {}).get('classification', "")
    return members, te_meta, raw_lines


# ============================================================
# Per-lineage clustering
# ============================================================

def get_seq(genome, chrom: str, start: int, end: int, strand: str) -> str:
    s = str(genome[chrom][start - 1:end]).upper()
    if strand == '-':
        s = revcomp(s)
    return s


def mmseqs_cluster_lineage(seqs: Dict[str, str], identity: float,
                           threads: int, tmpdir: str) -> Dict[str, str]:
    if len(seqs) < 2:
        return {k: k for k in seqs}
    fa = os.path.join(tmpdir, "in.fasta")
    with open(fa, "w") as f:
        for nm, s in seqs.items():
            f.write(f">{nm}\n{s}\n")
    pfx = os.path.join(tmpdir, "cl")
    mtmp = os.path.join(tmpdir, "mtmp")
    os.makedirs(mtmp, exist_ok=True)
    cmd = ["mmseqs", "easy-cluster", fa, pfx, mtmp,
           "--min-seq-id", str(identity),
           "-c", "0.8", "--cov-mode", "0",
           "--threads", str(threads), "-v", "0"]
    rc = subprocess.run(cmd, stdout=subprocess.DEVNULL,
                        stderr=subprocess.DEVNULL).returncode
    tsv = pfx + "_cluster.tsv"
    if rc != 0 or not os.path.exists(tsv) or os.path.getsize(tsv) == 0:
        return {k: k for k in seqs}
    out: Dict[str, str] = {}
    with open(tsv) as f:
        for line in f:
            rep, mem = line.rstrip().split('\t')
            out[mem] = rep
    return out


def build_clusters(members: List[Member], genome, identity: float,
                   threads: int) -> Tuple[List[List[Member]], List[str], List[str]]:
    """Cluster 5'LTR bodies per lineage via MMseqs2; for each cluster,
    add sibling 3'LTRs by parent.

    Returns (cluster_members, cluster_lineage, cluster_rep) where
    cluster_rep is the MMseqs2 representative member name (used as cluster_id).
    """
    by_parent_3: Dict[str, Member] = {}
    for m in members:
        if m.role == '3LTR':
            by_parent_3[m.parent] = m
    five_ltrs = [m for m in members if m.role == '5LTR']

    by_lineage: Dict[str, List[Member]] = defaultdict(list)
    for m in five_ltrs:
        if m.lineage:
            by_lineage[m.lineage].append(m)

    cluster_members: List[List[Member]] = []
    cluster_lineage: List[str] = []
    cluster_rep: List[str] = []
    tmp_root = tempfile.mkdtemp(prefix="dante_ltr_refine_clu_")
    try:
        for lineage, mems in by_lineage.items():
            seqs: Dict[str, str] = {}
            name_to_m: Dict[str, Member] = {}
            for m in mems:
                nm = f"{m.chrom}_{m.start}_{m.end}"
                seqs[nm] = get_seq(genome, m.chrom, m.start, m.end, m.strand)
                name_to_m[nm] = m
            ltmp = os.path.join(tmp_root,
                                re.sub(r"[^A-Za-z0-9]", "_", lineage))
            os.makedirs(ltmp, exist_ok=True)
            mship = mmseqs_cluster_lineage(seqs, identity, threads, ltmp)
            byrep: Dict[str, List[Member]] = defaultdict(list)
            for nm, rep in mship.items():
                if nm in name_to_m:
                    byrep[rep].append(name_to_m[nm])
            for rep, mem5 in byrep.items():
                joint: List[Member] = []
                for m in mem5:
                    joint.append(m)
                    sib = by_parent_3.get(m.parent)
                    if sib:
                        joint.append(sib)
                cluster_members.append(joint)
                cluster_lineage.append(lineage)
                cluster_rep.append(rep)
    finally:
        shutil.rmtree(tmp_root, ignore_errors=True)
    return cluster_members, cluster_lineage, cluster_rep


# ============================================================
# Per-element refinement record
# ============================================================

@dataclass
class RefRecord:
    """Per-LTR-feature refinement record (v2 schema).

    Coordinates are 1-based inclusive (GFF3 convention).  `start_orig`
    / `end_orig` are the DANTE_LTR input coords; `final_corrected_g`
    is the final refined OUTER-edge coordinate after policy + gates,
    or None when unrefined.
    """
    chrom: str
    start_orig: int
    end_orig: int
    strand: str
    role: str                  # '5LTR' | '3LTR'
    parent: str                # parent transposable_element ID
    lineage: str
    cluster_id: str = ""
    cluster_size: int = 0
    refinement_side: str = ""  # '5' (5LTR's outer 5') | '3' (3LTR's outer 3')

    # ---- per-pool diagnostics (raw, before gate) ----
    outer_corrected_g: Optional[int] = None
    outer_motif_ok: Optional[bool] = None
    outer_snap_offset: int = 0
    outer_n_pairs: int = 0

    inner_corrected_g: Optional[int] = None
    inner_motif_ok: Optional[bool] = None
    inner_snap_offset: int = 0
    inner_n_pairs: int = 0

    # ---- captured original state ----
    orig_motif_ok: Optional[bool] = None

    # ---- MAFFT fallback diagnostics (filled only if invoked) ----
    mafft_corrected_g: Optional[int] = None
    mafft_motif_ok: Optional[bool] = None

    # ---- final, after gates ----
    final_corrected_g: Optional[int] = None
    refinement_method: str = "none"   # parasail_inner | parasail_outer | mafft | none
    confidence: str = "unrefined"     # dual | divergent | inner_only | unrefined
    revert_reason: str = ""           # "" | "tsd_lost" | "tsd_lost_partner" | ...

    # ---- per-element TSD outcome (set on BOTH 5LTR and 3LTR records of an element) ----
    tsd_outcome: str = ""             # kept_exact | kept_fuzzy | shifted | lost | gained | both_none
    tsd_new_len: int = 0
    tsd_new_seq: str = ""
    orig_tsd_len: int = 0
    orig_tsd_seq: str = ""


def _classify_confidence(inner_g: Optional[int], inner_motif_ok: Optional[bool],
                         outer_g: Optional[int], outer_motif_ok: Optional[bool]
                         ) -> str:
    """Confidence label policy (v2).

    dual       : inner accepted, outer also reached a valid coord, |outer-inner| <= 5
    divergent  : inner accepted, outer also reached a valid coord, |outer-inner| > 5
    inner_only : inner accepted, outer did not validate
    unrefined  : inner did not validate (we keep the original)

    Caller has already decided to use the inner coord (this is invoked
    only when inner_motif_ok is True); for `unrefined` (inner failed),
    pass inner_motif_ok=False/None and we return "unrefined".
    """
    if inner_motif_ok is not True or inner_g is None:
        return "unrefined"
    if outer_motif_ok is True and outer_g is not None:
        return "dual" if abs(outer_g - inner_g) <= 5 else "divergent"
    return "inner_only"


# ============================================================
# Cluster orchestration
# ============================================================

def _outer_g(member: Member, side: str) -> int:
    """Genomic coord of the outer boundary for `side` on `member`."""
    if (side == "5" and member.strand == "+") or (side == "3" and member.strand == "-"):
        return member.start
    return member.end


def refine_one_cluster(cluster_members: List[Member], cluster_id: str,
                       cluster_lineage: str,
                       genome, genome_lens: Dict[str, int],
                       opts: argparse.Namespace
                       ) -> List[RefRecord]:
    """v2: run two parasail passes per side (outer-edge + inner-edge pool),
    apply the inner-primary policy, assign confidence labels.

    Per-element TSD recheck and revert is applied later in `refine_all`
    after MAFFT fallback (so the gate sees the final coord set).
    Returns one RefRecord per LTR feature in the cluster.
    """
    members_5 = [m for m in cluster_members if m.role == "5LTR"]
    members_3 = [m for m in cluster_members if m.role == "3LTR"]
    if not members_5 or not members_3:
        # Cluster lacks one role entirely; cannot run inner-pool.
        return [
            RefRecord(chrom=m.chrom, start_orig=m.start, end_orig=m.end,
                      strand=m.strand, role=m.role, parent=m.parent,
                      lineage=m.lineage, cluster_id=cluster_id,
                      cluster_size=len(cluster_members),
                      refinement_side="5" if m.role == "5LTR" else "3")
            for m in cluster_members
        ]

    # Pre-extract anchor pools per side.  extract_seq_for_side() uses
    # coordinates + strand only (role label irrelevant), so calling
    # it on opposite-role members yields the inner-edge geometry we
    # want.  See docs/refine_v2_analysis.md §2.
    pool_outer_5 = extract_pool_for_side(members_5, "5LTR", "5", genome,
                                         genome_lens, opts.anchor_len, opts.flank_len)
    pool_inner_5 = extract_pool_for_side(members_3, "3LTR", "5", genome,
                                         genome_lens, opts.anchor_len, opts.flank_len)
    pool_outer_3 = extract_pool_for_side(members_3, "3LTR", "3", genome,
                                         genome_lens, opts.anchor_len, opts.flank_len)
    pool_inner_3 = extract_pool_for_side(members_5, "5LTR", "3", genome,
                                         genome_lens, opts.anchor_len, opts.flank_len)

    # Index queries by (member feat_id, side) so we can look up by member.
    queries_5 = {q["member"].feat_id: q for q in pool_outer_5}
    queries_3 = {q["member"].feat_id: q for q in pool_outer_3}

    records: List[RefRecord] = []
    for m in cluster_members:
        rec = RefRecord(
            chrom=m.chrom, start_orig=m.start, end_orig=m.end,
            strand=m.strand, role=m.role, parent=m.parent, lineage=m.lineage,
            cluster_id=cluster_id,
            cluster_size=len(cluster_members),
            refinement_side="5" if m.role == "5LTR" else "3",
        )

        side = rec.refinement_side
        # Capture original motif at the original outer-boundary coord.
        chrom_len = genome_lens[m.chrom]
        og = _outer_g(m, side)
        orig_pair = get_boundary_pair(genome, m.chrom, m.strand, og, side, chrom_len)
        rec.orig_motif_ok = validate_motif(orig_pair, side, "TG/CA")

        if side == "5":
            qmap = queries_5
            outer_pool = [a for a in pool_outer_5 if a["member"].feat_id != m.feat_id]
            inner_pool = pool_inner_5  # cannot include self (opposite-role)
        else:
            qmap = queries_3
            outer_pool = [a for a in pool_outer_3 if a["member"].feat_id != m.feat_id]
            inner_pool = pool_inner_3
        q = qmap.get(m.feat_id)
        if q is None:
            records.append(rec)
            continue

        # Outer-pool extension
        outer_res = aggregate_extension(
            q, outer_pool, side,
            match=opts.match, mismatch=opts.mismatch,
            gap_open=opts.gap_open, gap_extend=opts.gap_extend,
            score_threshold=opts.score_threshold,
            min_num_alignments=opts.min_num_alignments,
            max_pairs=opts.max_pairs,
        )
        outer_g = project_corrected_g(m, side, outer_res["ext_len"], q["anchor"])
        if outer_g is not None:
            outer_g_snap, outer_off, outer_ok = snap_to_motif(
                genome, m, outer_g, side, chrom_len, opts.snap_window, "TG/CA")
        else:
            outer_g_snap, outer_off, outer_ok = (None, 0, None)
        rec.outer_corrected_g = outer_g_snap
        rec.outer_motif_ok = outer_ok
        rec.outer_snap_offset = outer_off
        rec.outer_n_pairs = outer_res["n_pairs"]

        # Inner-pool extension
        inner_res = aggregate_extension(
            q, inner_pool, side,
            match=opts.match, mismatch=opts.mismatch,
            gap_open=opts.gap_open, gap_extend=opts.gap_extend,
            score_threshold=opts.score_threshold,
            min_num_alignments=opts.min_num_alignments,
            max_pairs=opts.max_pairs,
        )
        inner_g = project_corrected_g(m, side, inner_res["ext_len"], q["anchor"])
        if inner_g is not None:
            inner_g_snap, inner_off, inner_ok = snap_to_motif(
                genome, m, inner_g, side, chrom_len, opts.snap_window, "TG/CA")
        else:
            inner_g_snap, inner_off, inner_ok = (None, 0, None)
        rec.inner_corrected_g = inner_g_snap
        rec.inner_motif_ok = inner_ok
        rec.inner_snap_offset = inner_off
        rec.inner_n_pairs = inner_res["n_pairs"]

        # Apply inner-primary policy.
        if inner_ok is True and inner_g_snap is not None:
            rec.final_corrected_g = inner_g_snap
            rec.refinement_method = "parasail_inner"
            rec.confidence = _classify_confidence(
                inner_g_snap, inner_ok, outer_g_snap, outer_ok)
        # else: inner failed; leave final_corrected_g=None (unrefined).

        records.append(rec)
    return records


def cluster_inner_validation_rate(records: List[RefRecord]
                                   ) -> Tuple[float, int, int]:
    """Fraction of records in this cluster where the INNER pool validated.

    Returns (rate, n_validated, n_total).
    """
    n_total = len(records)
    if n_total == 0:
        return 1.0, 0, 0
    n_validated = sum(1 for r in records if r.inner_motif_ok is True)
    return n_validated / n_total, n_validated, n_total


def cluster_outer_validation_rate(records: List[RefRecord]
                                   ) -> Tuple[float, int, int]:
    """Fraction of records in this cluster where the OUTER pool validated."""
    n_total = len(records)
    if n_total == 0:
        return 1.0, 0, 0
    n_validated = sum(1 for r in records if r.outer_motif_ok is True)
    return n_validated / n_total, n_validated, n_total


def apply_tsd_gate_per_element(records: List[RefRecord], te_meta: Dict[str, Dict],
                                 genome, genome_lens: Dict[str, int],
                                 enable_revert: bool = True) -> None:
    """Per-element TSD recheck + revert rule.

    For each TE with both LTRs in `records`, compute the TSD at the
    proposed coordinates.  If the original had a TSD AND the new
    coordinates produce a `lost` outcome, revert all *changed* sides
    on that element to DANTE_LTR original (set final_corrected_g back
    to None, refinement_method='none', confidence='unrefined',
    revert_reason='tsd_lost').

    Always populates rec.tsd_outcome / tsd_new_len / tsd_new_seq /
    orig_tsd_len / orig_tsd_seq on each member of an element.

    `enable_revert=False` keeps the proposed coords (diagnostic mode);
    the tsd_outcome attribute is still populated.
    """
    by_parent: Dict[str, Dict[str, RefRecord]] = defaultdict(dict)
    for r in records:
        by_parent[r.parent][r.role] = r

    for te_id, pair in by_parent.items():
        rec5 = pair.get("5LTR"); rec3 = pair.get("3LTR")
        meta = te_meta.get(te_id, {})
        orig_tsd_seq = meta.get("tsd_seq", "")
        orig_tsd_len = int(meta.get("tsd_len", 0) or 0)
        orig_had_tsd = orig_tsd_len >= 4

        # Proposed coords (fall back to original where side wasn't refined)
        if rec5 is not None:
            v5 = rec5.final_corrected_g if rec5.final_corrected_g is not None \
                 else _outer_g(_member_from_rec(rec5), "5")
        if rec3 is not None:
            v3 = rec3.final_corrected_g if rec3.final_corrected_g is not None \
                 else _outer_g(_member_from_rec(rec3), "3")
        if rec5 is None or rec3 is None:
            # Can't TSD-check without both sides — skip; populate orig only.
            for r in (rec5, rec3):
                if r is None:
                    continue
                r.orig_tsd_seq = orig_tsd_seq
                r.orig_tsd_len = orig_tsd_len
            continue

        chrom = rec5.chrom
        strand = rec5.strand
        new_n, new_seq, fuzzy = detect_tsd(genome, chrom, v5, v3, strand)
        new_has = new_n >= 4

        # Categorize outcome
        if not orig_had_tsd and not new_has:
            outcome = "both_none"
        elif orig_had_tsd and new_has:
            o = orig_tsd_seq.upper()
            n = new_seq.upper()
            if "/" in o:
                halves = o.split("/", 1)
                outcome = "kept_fuzzy" if n in halves else "shifted"
            else:
                outcome = "kept_exact" if n == o else "shifted"
        elif orig_had_tsd and not new_has:
            outcome = "lost"
        else:
            outcome = "gained"

        for r in (rec5, rec3):
            r.tsd_outcome = outcome
            r.tsd_new_len = new_n
            r.tsd_new_seq = new_seq + (";fuzzy" if fuzzy else "")
            r.orig_tsd_seq = orig_tsd_seq
            r.orig_tsd_len = orig_tsd_len

        # Apply revert rule: lost AND original had TSD → revert any changed side
        if enable_revert and outcome == "lost" and orig_had_tsd:
            for r in (rec5, rec3):
                if r.final_corrected_g is not None:
                    r.final_corrected_g = None
                    r.refinement_method = "none"
                    r.confidence = "unrefined"
                    r.revert_reason = "tsd_lost"
            # Re-detect TSD at original coords for reporting
            v5_orig = _outer_g(_member_from_rec(rec5), "5")
            v3_orig = _outer_g(_member_from_rec(rec3), "3")
            on, oseq, ofz = detect_tsd(genome, chrom, v5_orig, v3_orig, strand)
            for r in (rec5, rec3):
                r.tsd_outcome = "kept_exact" if on >= 4 else "both_none"
                r.tsd_new_len = on
                r.tsd_new_seq = oseq + (";fuzzy" if ofz else "")


def _member_from_rec(rec: RefRecord) -> Member:
    """Build a Member from a RefRecord (for outer-coord helpers)."""
    return Member(feat_id=0, chrom=rec.chrom,
                  start=rec.start_orig, end=rec.end_orig,
                  strand=rec.strand, role=rec.role,
                  parent=rec.parent, lineage=rec.lineage)


# ============================================================
# MAFFT fallback (subprocess to refine_mafft_fallback.R)
# ============================================================

def run_mafft_fallback(records: List[RefRecord], cluster_members: List[Member],
                       cluster_id: str, cluster_lineage: str,
                       genome, genome_lens: Dict[str, int],
                       opts: argparse.Namespace) -> None:
    """Mutate `records` in place with MAFFT-corrected boundaries for any
    record whose parasail call did NOT validate, IFF MAFFT's per-member
    call validates better.

    Uses an Rscript subprocess.  Writes a per-cluster FASTA of LTR ±
    --mafft_flank bp and reads back a per-member TSV.
    """
    script = os.path.join(HERE, "refine_mafft_fallback.R")
    if not os.path.exists(script):
        logger.warning("MAFFT fallback script missing at %s; skipping", script)
        return

    tmpdir = tempfile.mkdtemp(prefix="dante_ltr_refine_mafft_")
    try:
        fa_path = os.path.join(tmpdir, "cluster.fasta")
        meta_path = os.path.join(tmpdir, "members.tsv")
        out_path = os.path.join(tmpdir, "result.tsv")

        # Build per-member extended FASTA in biological orientation.
        # Each row is: header = "<chrom>:<start>-<end>:<role>:<strand>"
        # (so we can match back).
        with open(fa_path, "w") as fa, open(meta_path, "w") as meta:
            meta.write("name\tchrom\tstart\tend\tstrand\trole\t"
                       "ext_start\text_end\tbody_start_raw\tbody_end_raw\n")
            for m in cluster_members:
                glen = genome_lens[m.chrom]
                ext_start = max(1, m.start - opts.mafft_flank)
                ext_end = min(glen, m.end + opts.mafft_flank)
                seq = str(genome[m.chrom][ext_start - 1:ext_end]).upper()
                if m.strand == '-':
                    seq = revcomp(seq)
                # In bio orientation, body offsets:
                gleft = m.start - ext_start
                gright = ext_end - m.end
                bio_5flank = gright if m.strand == '-' else gleft
                bio_body = m.end - m.start + 1
                body_start_raw = bio_5flank + 1
                body_end_raw = bio_5flank + bio_body
                name = f"{m.chrom}__{m.start}__{m.end}__{m.role}__{m.strand}"
                fa.write(f">{name}\n{seq}\n")
                meta.write(f"{name}\t{m.chrom}\t{m.start}\t{m.end}\t"
                           f"{m.strand}\t{m.role}\t{ext_start}\t{ext_end}\t"
                           f"{body_start_raw}\t{body_end_raw}\n")

        cmd = ["Rscript", script,
               "--input_fasta", fa_path,
               "--members_tsv", meta_path,
               "--output_tsv", out_path,
               "--flank", str(opts.mafft_flank),
               "--threads", str(opts.threads),
               "--boundary_motif", "TG/CA"]
        try:
            r = subprocess.run(cmd, capture_output=True, text=True,
                               check=False, timeout=opts.mafft_timeout)
        except subprocess.TimeoutExpired:
            logger.warning("MAFFT fallback timeout for cluster %s", cluster_id)
            return
        if r.returncode != 0:
            logger.warning("MAFFT fallback failed for cluster %s: %s",
                           cluster_id, r.stderr.strip()[:200])
            return
        if not os.path.exists(out_path):
            return

        # Parse: per-row member name + corrected_g + motif_ok
        by_name: Dict[str, Tuple[Optional[int], Optional[bool]]] = {}
        with open(out_path) as f:
            header = f.readline().rstrip().split("\t")
            try:
                i_name = header.index("name")
                i_g = header.index("corrected_g")
                i_ok = header.index("motif_ok")
            except ValueError:
                return
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                nm = parts[i_name]
                cg = parts[i_g]
                mo = parts[i_ok]
                cg_v = int(cg) if cg not in ("NA", "", "None") else None
                mo_v = (True if mo == "TRUE" else
                        False if mo == "FALSE" else None)
                by_name[nm] = (cg_v, mo_v)

        # Merge into records.  Override only on records the inner-primary
        # policy did NOT already accept (final_corrected_g is None) AND
        # where MAFFT's call validates the TG/CA motif.  TSD-loss revert
        # is applied later in refine_all by `apply_tsd_gate_per_element`.
        for rec in records:
            nm = (f"{rec.chrom}__{rec.start_orig}__{rec.end_orig}"
                  f"__{rec.role}__{rec.strand}")
            mafft = by_name.get(nm)
            if mafft is None:
                continue
            cg, mok = mafft
            rec.mafft_corrected_g = cg
            rec.mafft_motif_ok = mok
            if rec.final_corrected_g is not None:
                continue   # parasail_inner already accepted this side
            if cg is None or mok is not True:
                continue
            rec.final_corrected_g = cg
            rec.refinement_method = "mafft"
            rec.confidence = "inner_only"
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


# ============================================================
# GFF3 emission
# ============================================================

def _strip_attrs(attrs: str, keys_to_drop: Tuple[str, ...]) -> str:
    """Remove keys from a GFF3 attribute string, preserving order."""
    kept = []
    for k, v in _ATTR_RE.findall(attrs):
        if k in keys_to_drop:
            continue
        kept.append(f"{k}={v}")
    return ";".join(kept)


def _refined_ltr_coords(rec: RefRecord) -> Tuple[int, int, str]:
    """Return (new_start, new_end, applied) for an LTR feature given its
    refinement record.

    `applied` is "yes" when the start/end actually changed and
    refinement_method != "none", else "no".
    """
    s, e = rec.start_orig, rec.end_orig
    if rec.refinement_method == "none" or rec.final_corrected_g is None:
        return s, e, "no"

    cg = rec.final_corrected_g
    if rec.role == "5LTR":
        if rec.strand == '+':
            new_s, new_e = cg, e
        else:
            new_s, new_e = s, cg
    else:  # 3LTR
        if rec.strand == '+':
            new_s, new_e = s, cg
        else:
            new_s, new_e = cg, e

    # Sanity guard: keep start <= end, fall back if not
    if new_s > new_e:
        return s, e, "no"
    return new_s, new_e, "yes"


_CONF_RANK = {"unrefined": 0, "inner_only": 1, "divergent": 2, "dual": 3}


def _tsd_child_coords(rec5: RefRecord, rec3: RefRecord, tsd_len: int
                      ) -> Tuple[Tuple[int, int], Tuple[int, int]]:
    """Genomic coords of (TSD_L, TSD_R) child rows given final outer
    boundaries and a TSD length.  Returns ((l_start,l_end),(r_start,r_end))
    in 1-based inclusive forward strand convention.
    """
    v5 = (rec5.final_corrected_g if rec5.final_corrected_g is not None
          else _outer_g(_member_from_rec(rec5), "5"))
    v3 = (rec3.final_corrected_g if rec3.final_corrected_g is not None
          else _outer_g(_member_from_rec(rec3), "3"))
    if rec5.strand == "+":
        l_start, l_end = v5 - tsd_len, v5 - 1
        r_start, r_end = v3 + 1, v3 + tsd_len
    else:
        # On - strand, biological upstream of element = + strand right of v5;
        # biological downstream = + strand left of v3.
        l_start, l_end = v5 + 1, v5 + tsd_len
        r_start, r_end = v3 - tsd_len, v3 - 1
    return (l_start, l_end), (r_start, r_end)


def emit_refined_gff3(raw_lines: List[str], records: List[RefRecord],
                      out_path: str) -> Dict[str, int]:
    """Emit refined GFF3 by replaying input lines and updating
    TE / LTR / TSD rows.

    For each transposable_element:
      - the TE row's coords are envelope of refined LTR coords;
      - the LTR rows have refined outer coords;
      - input target_site_duplication rows are dropped; new TSD rows
        are emitted immediately after the TE row, based on the
        per-element tsd_outcome.

    Returns counts dict for logging.
    """
    rec_idx: Dict[Tuple[str, int, int, str], RefRecord] = {}
    for r in records:
        rec_idx[(r.chrom, r.start_orig, r.end_orig, r.role)] = r
    rec_by_parent: Dict[str, Dict[str, RefRecord]] = defaultdict(dict)
    for r in records:
        rec_by_parent[r.parent][r.role] = r

    counts = {"te_total": 0, "te_refined": 0,
              "ltr_total": 0, "ltr_refined": 0,
              "tsd_emitted": 0, "tsd_dropped": 0}

    drop_keys = ("Original_Start", "Original_End",
                 "Refinement_Method", "Refinement_Confidence",
                 "Cluster_ID", "Cluster_Size", "TG_OK", "CA_OK",
                 "Motif_Orig", "Motif_New", "TSD_Outcome",
                 "Outer_Pool_g", "Outer_Pool_Motif_OK",
                 "Inner_Pool_g", "Inner_Pool_Motif_OK",
                 "Revert_Reason", "TSD_Status")

    with open(out_path, "w") as out:
        for line in raw_lines:
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                out.write(line)
                continue
            parts = stripped.split('\t')
            if len(parts) < 9:
                out.write(line)
                continue
            seqname, src, ftype, s1, e1, sc, strand, ph, attrs = parts
            attrs_clean = _strip_attrs(attrs, drop_keys)

            if ftype == 'target_site_duplication':
                # Drop input TSD rows; we emit fresh ones after the TE row.
                counts["tsd_dropped"] += 1
                continue

            if ftype == 'long_terminal_repeat':
                counts["ltr_total"] += 1
                role = _attr(attrs, 'LTR') or ''
                key = (seqname, int(s1), int(e1), role)
                rec = rec_idx.get(key)
                if rec is None:
                    out.write(line)
                    continue
                ns, ne, applied = _refined_ltr_coords(rec)
                if applied == "yes":
                    counts["ltr_refined"] += 1
                # Per-side motif at orig and new coords (2-bp dinucleotide).
                motif_orig = "TG" if (role == "5LTR" and rec.orig_motif_ok is True) \
                              else ("CA" if (role == "3LTR" and rec.orig_motif_ok is True)
                                    else "none")
                if rec.final_corrected_g is not None:
                    target = "TG" if role == "5LTR" else "CA"
                    motif_new = target  # accepted by inner-primary policy => motif_ok
                elif rec.orig_motif_ok is True:
                    motif_new = motif_orig
                else:
                    motif_new = "none"
                add_attrs_parts = [
                    f"Refinement_Method={rec.refinement_method}",
                    f"Refinement_Confidence={rec.confidence}",
                    f"Cluster_ID={rec.cluster_id}",
                    f"Cluster_Size={rec.cluster_size}",
                    f"Motif_Orig={motif_orig}",
                    f"Motif_New={motif_new}",
                ]
                if applied == "yes" or (ns != int(s1) or ne != int(e1)):
                    add_attrs_parts.insert(0, f"Original_End={rec.end_orig}")
                    add_attrs_parts.insert(0, f"Original_Start={rec.start_orig}")
                if rec.outer_corrected_g is not None:
                    add_attrs_parts.append(
                        f"Outer_Pool_g={rec.outer_corrected_g};"
                        f"Outer_Pool_Motif_OK={_b(rec.outer_motif_ok)}")
                if rec.inner_corrected_g is not None:
                    add_attrs_parts.append(
                        f"Inner_Pool_g={rec.inner_corrected_g};"
                        f"Inner_Pool_Motif_OK={_b(rec.inner_motif_ok)}")
                if rec.revert_reason:
                    add_attrs_parts.append(f"Revert_Reason={rec.revert_reason}")
                add_attrs = ";".join(add_attrs_parts)
                new_attrs = attrs_clean + ";" + add_attrs if attrs_clean else add_attrs
                out.write("\t".join([seqname, src, ftype, str(ns), str(ne),
                                     sc, strand, ph, new_attrs]) + "\n")

            elif ftype == 'transposable_element':
                counts["te_total"] += 1
                te_id = _attr(attrs, 'ID') or ''
                pair = rec_by_parent.get(te_id, {})
                rec5 = pair.get("5LTR"); rec3 = pair.get("3LTR")

                # Envelope from refined LTR coords; fall back to orig.
                refined_any = False
                env_s, env_e = int(s1), int(e1)
                if rec5 is not None or rec3 is not None:
                    ns5, ne5, ap5 = (_refined_ltr_coords(rec5) if rec5 is not None
                                      else (int(s1), int(e1), "no"))
                    ns3, ne3, ap3 = (_refined_ltr_coords(rec3) if rec3 is not None
                                      else (int(s1), int(e1), "no"))
                    if ap5 == "yes" or ap3 == "yes":
                        refined_any = True
                    if rec5 is not None and rec3 is not None:
                        env_s = min(ns5, ns3); env_e = max(ne5, ne3)
                    else:
                        env_s = min(ns5, ns3, int(s1))
                        env_e = max(ne5, ne3, int(e1))

                if refined_any:
                    counts["te_refined"] += 1

                # Element-level method/confidence: take the lower of the two sides.
                m5 = rec5.refinement_method if rec5 else "none"
                m3 = rec3.refinement_method if rec3 else "none"
                te_method = (m5 if m5 != "none" else m3 if m3 != "none" else "none")
                c5 = rec5.confidence if rec5 else "unrefined"
                c3 = rec3.confidence if rec3 else "unrefined"
                te_conf = min((c5, c3), key=lambda c: _CONF_RANK.get(c, 0))
                cluster_id = (rec5.cluster_id if rec5 else
                              rec3.cluster_id if rec3 else "")
                cluster_size = (rec5.cluster_size if rec5 else
                                rec3.cluster_size if rec3 else 0)

                # TSD outcome (set on either rec5 or rec3 — both carry the same value)
                tsd_outcome = ""
                tsd_new_len = 0
                tsd_new_seq = ""
                for r in (rec5, rec3):
                    if r is not None and r.tsd_outcome:
                        tsd_outcome = r.tsd_outcome
                        tsd_new_len = r.tsd_new_len
                        tsd_new_seq = r.tsd_new_seq
                        break

                # Decide TE-row TSD= attribute value
                if tsd_outcome in ("kept_exact", "kept_fuzzy", "shifted", "gained") \
                        and tsd_new_len >= 4:
                    tsd_attr_val = (tsd_new_seq.split(";", 1)[0]
                                    if ";fuzzy" in tsd_new_seq else tsd_new_seq)
                elif tsd_outcome == "lost":
                    tsd_attr_val = "not_found"
                else:
                    # both_none or no element-pair info
                    tsd_attr_val = (tsd_new_seq if tsd_new_len >= 4 else "not_found")

                # Drop input TSD= attr; we'll re-add ours.
                attrs_clean_te = _strip_attrs(attrs_clean, ("TSD",))
                # Per-side motif summary (5LTR motif = TG; 3LTR = CA)
                motif_orig_5 = "TG" if (rec5 and rec5.orig_motif_ok is True) else "none"
                motif_new_5 = "TG" if (rec5 and (rec5.final_corrected_g is not None
                                                  or rec5.orig_motif_ok is True)) else "none"
                motif_orig_3 = "CA" if (rec3 and rec3.orig_motif_ok is True) else "none"
                motif_new_3 = "CA" if (rec3 and (rec3.final_corrected_g is not None
                                                  or rec3.orig_motif_ok is True)) else "none"

                add_parts = [
                    f"Original_Start={s1}",
                    f"Original_End={e1}",
                    f"TSD={tsd_attr_val}",
                    f"Refinement_Method={te_method}",
                    f"Refinement_Confidence={te_conf}",
                    f"Cluster_ID={cluster_id}",
                    f"Cluster_Size={cluster_size}",
                    f"Motif_Orig={motif_orig_5}/{motif_orig_3}",
                    f"Motif_New={motif_new_5}/{motif_new_3}",
                ]
                if tsd_outcome:
                    add_parts.append(f"TSD_Outcome={tsd_outcome}")
                if rec5 is not None and rec5.outer_corrected_g is not None \
                        and rec3 is not None and rec3.outer_corrected_g is not None:
                    add_parts.append(
                        f"Outer_Pool_g={rec5.outer_corrected_g},{rec3.outer_corrected_g}")
                    add_parts.append(
                        f"Outer_Pool_Motif_OK={_b(rec5.outer_motif_ok)},"
                        f"{_b(rec3.outer_motif_ok)}")
                if rec5 is not None and rec5.inner_corrected_g is not None \
                        and rec3 is not None and rec3.inner_corrected_g is not None:
                    add_parts.append(
                        f"Inner_Pool_g={rec5.inner_corrected_g},{rec3.inner_corrected_g}")
                    add_parts.append(
                        f"Inner_Pool_Motif_OK={_b(rec5.inner_motif_ok)},"
                        f"{_b(rec3.inner_motif_ok)}")

                add_attrs = ";".join(add_parts)
                new_attrs = attrs_clean_te + ";" + add_attrs if attrs_clean_te else add_attrs
                out.write("\t".join([seqname, src, ftype, str(env_s), str(env_e),
                                     sc, strand, ph, new_attrs]) + "\n")

                # Emit fresh target_site_duplication rows immediately after the TE.
                if rec5 is not None and rec3 is not None and tsd_new_len >= 4 \
                        and tsd_outcome in ("kept_exact", "kept_fuzzy", "shifted", "gained"):
                    (l_s, l_e), (r_s, r_e) = _tsd_child_coords(rec5, rec3, tsd_new_len)
                    tsd_status = ("original" if tsd_outcome == "kept_exact" else
                                  "re_identified" if tsd_outcome in ("kept_fuzzy", "shifted")
                                  else "gained")
                    common = (f"Parent={te_id};TSD_Status={tsd_status};"
                              f"Refinement_Confidence={te_conf}")
                    out.write("\t".join([
                        seqname, src, "target_site_duplication",
                        str(l_s), str(l_e), ".", strand, ".", common]) + "\n")
                    out.write("\t".join([
                        seqname, src, "target_site_duplication",
                        str(r_s), str(r_e), ".", strand, ".", common]) + "\n")
                    counts["tsd_emitted"] += 2
            else:
                out.write(line)
    return counts


def _b(v: Optional[bool]) -> str:
    if v is True: return "TRUE"
    if v is False: return "FALSE"
    return "NA"


# ============================================================
# Per-element TSV
# ============================================================

def emit_per_element_tsv(records: List[RefRecord], out_path: str) -> None:
    cols = [
        "chrom", "start_orig", "end_orig", "strand", "role",
        "parent_te_id", "lineage_full",
        "cluster_id", "cluster_size", "refinement_side",
        # outer pool
        "outer_corrected_g", "outer_motif_ok", "outer_snap_offset", "outer_n_pairs",
        # inner pool
        "inner_corrected_g", "inner_motif_ok", "inner_snap_offset", "inner_n_pairs",
        # captured original state
        "orig_motif_ok",
        # mafft fallback (filled only when invoked)
        "mafft_corrected_g", "mafft_motif_ok",
        # final
        "final_corrected_g", "refinement_method", "confidence", "revert_reason",
        # tsd outcome
        "orig_tsd_len", "orig_tsd_seq",
        "tsd_outcome", "tsd_new_len", "tsd_new_seq",
        # convenience
        "shift_bp",
    ]

    def fmt(v):
        if v is None:
            return "NA"
        if isinstance(v, bool):
            return "TRUE" if v else "FALSE"
        return str(v)

    with open(out_path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in records:
            anchored = _outer_g(_member_from_rec(r), r.refinement_side)
            shift = (r.final_corrected_g - anchored
                     if r.final_corrected_g is not None else None)
            row = [
                r.chrom, r.start_orig, r.end_orig, r.strand, r.role,
                r.parent, r.lineage,
                r.cluster_id, r.cluster_size, r.refinement_side,
                r.outer_corrected_g, r.outer_motif_ok, r.outer_snap_offset, r.outer_n_pairs,
                r.inner_corrected_g, r.inner_motif_ok, r.inner_snap_offset, r.inner_n_pairs,
                r.orig_motif_ok,
                r.mafft_corrected_g, r.mafft_motif_ok,
                r.final_corrected_g, r.refinement_method, r.confidence, r.revert_reason,
                r.orig_tsd_len, r.orig_tsd_seq,
                r.tsd_outcome, r.tsd_new_len, r.tsd_new_seq,
                shift,
            ]
            f.write("\t".join(fmt(v) for v in row) + "\n")


# ============================================================
# Cluster manifest
# ============================================================

def emit_cluster_manifest(records: List[RefRecord], cluster_index,
                          out_path: str) -> None:
    """Cluster manifest TSV.

    cluster_index: dict cluster_id -> dict with 'lineage', 'mafft_invoked',
    'mafft_validation_rate'.
    """
    by_cl: Dict[str, List[RefRecord]] = defaultdict(list)
    for r in records:
        if r.cluster_id:
            by_cl[r.cluster_id].append(r)

    with open(out_path, "w") as f:
        f.write("cluster_id\tlineage_full\tn_total\tn_5ltr\tn_3ltr\t"
                "inner_validated_5\tinner_validated_3\t"
                "outer_validated_5\touter_validated_3\t"
                "inner_validation_rate\touter_validation_rate\t"
                "mafft_invoked\tmafft_validation_rate\tlow_confidence\n")
        for cid, recs in by_cl.items():
            lin = recs[0].lineage if recs else ""
            n5 = sum(1 for r in recs if r.role == "5LTR")
            n3 = sum(1 for r in recs if r.role == "3LTR")
            iv5 = sum(1 for r in recs if r.role == "5LTR" and r.inner_motif_ok is True)
            iv3 = sum(1 for r in recs if r.role == "3LTR" and r.inner_motif_ok is True)
            ov5 = sum(1 for r in recs if r.role == "5LTR" and r.outer_motif_ok is True)
            ov3 = sum(1 for r in recs if r.role == "3LTR" and r.outer_motif_ok is True)
            denom = len(recs) or 1
            inner_rate = (iv5 + iv3) / denom
            outer_rate = (ov5 + ov3) / denom
            cmeta = cluster_index.get(cid, {})
            mafft_inv = "TRUE" if cmeta.get("mafft_invoked") else "FALSE"
            mafft_rate = cmeta.get("mafft_validation_rate")
            mafft_rate_str = f"{mafft_rate:.4f}" if mafft_rate is not None else "NA"
            low_conf = "TRUE" if (iv5 + iv3) < 4 else "FALSE"
            f.write(f"{cid}\t{lin}\t{len(recs)}\t{n5}\t{n3}\t"
                    f"{iv5}\t{iv3}\t{ov5}\t{ov3}\t"
                    f"{inner_rate:.4f}\t{outer_rate:.4f}\t"
                    f"{mafft_inv}\t{mafft_rate_str}\t{low_conf}\n")


# ============================================================
# Main driver
# ============================================================

def refine_all(args: argparse.Namespace) -> Dict:
    from pyfaidx import Fasta  # type: ignore

    out_prefix = args.output
    out_dir = os.path.dirname(out_prefix) or "."
    if out_dir != ".":
        os.makedirs(out_dir, exist_ok=True)

    t0_total = time.time()

    logger.info("Reading GFF3: %s", args.gff3)
    members, te_meta, raw_lines = parse_gff3_for_refinement(args.gff3)
    logger.info("  %d LTR features, %d TEs", len(members), len(te_meta))

    logger.info("Opening genome: %s", args.genome)
    genome = Fasta(args.genome, sequence_always_upper=True, as_raw=True)
    genome_lens = {n: len(genome[n]) for n in genome.keys()}

    logger.info("Clustering 5'LTR bodies per lineage at identity=%.2f",
                args.identity)
    t0 = time.time()
    cluster_members, cluster_lineage, cluster_rep = build_clusters(
        members, genome, args.identity, args.threads)
    t_cluster = time.time() - t0
    logger.info("  built %d clusters in %.1fs", len(cluster_members), t_cluster)

    # Per-role cluster-size filter (v2): require N(5LTR) >= N AND N(3LTR) >= N.
    def _qualifies(mems):
        n5 = sum(1 for m in mems if m.role == "5LTR")
        n3 = sum(1 for m in mems if m.role == "3LTR")
        return n5 >= args.min_cluster_size and n3 >= args.min_cluster_size
    qualifying = [(ci, mems, cluster_lineage[ci], cluster_rep[ci])
                  for ci, mems in enumerate(cluster_members)
                  if _qualifies(mems)]
    logger.info("  %d clusters qualify (>= --min_cluster_size=%d of EACH role)",
                len(qualifying), args.min_cluster_size)

    t0 = time.time()
    cluster_records: Dict[str, List[RefRecord]] = {}
    cluster_index: Dict[str, Dict] = {}

    def _do_one(task):
        ci, mems, lin, rep = task
        cluster_id = f"clu_{ci:06d}_{rep[:80]}"
        recs = refine_one_cluster(mems, cluster_id, lin,
                                  genome, genome_lens, args)
        return cluster_id, lin, recs

    if args.workers > 1 and len(qualifying) > 1:
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futs = [pool.submit(_do_one, t) for t in qualifying]
            done = 0
            for fut in as_completed(futs):
                cid, lin, recs = fut.result()
                cluster_records[cid] = recs
                cluster_index[cid] = {"lineage": lin, "mafft_invoked": False,
                                      "mafft_validation_rate": None}
                done += 1
                if done == 1 or done % 50 == 0:
                    logger.info("  parasail [%d/%d]", done, len(qualifying))
    else:
        for k, t in enumerate(qualifying, 1):
            cid, lin, recs = _do_one(t)
            cluster_records[cid] = recs
            cluster_index[cid] = {"lineage": lin, "mafft_invoked": False,
                                  "mafft_validation_rate": None}
            if k == 1 or k % 50 == 0:
                logger.info("  parasail [%d/%d]", k, len(qualifying))
    t_parasail = time.time() - t0
    logger.info("Parasail finished in %.1fs", t_parasail)

    # MAFFT fallback per cluster (gated on inner-pool validation rate).
    t0 = time.time()
    n_fallback = 0
    if args.mafft_fallback:
        cl_mems_by_id: Dict[str, List[Member]] = {}
        for ci, mems in enumerate(cluster_members):
            if not _qualifies(mems):
                continue
            cluster_id = f"clu_{ci:06d}_{cluster_rep[ci][:80]}"
            cl_mems_by_id[cluster_id] = mems
        for cid, recs in cluster_records.items():
            rate, _, _ = cluster_inner_validation_rate(recs)
            if rate < args.mafft_fallback_threshold:
                n_fallback += 1
                cluster_index[cid]["mafft_invoked"] = True
                run_mafft_fallback(recs, cl_mems_by_id[cid],
                                   cid, cluster_index[cid]["lineage"],
                                   genome, genome_lens, args)
                n_accepted = sum(1 for r in recs if r.final_corrected_g is not None)
                cluster_index[cid]["mafft_validation_rate"] = (
                    n_accepted / len(recs) if recs else 0.0)
    t_mafft = time.time() - t0
    logger.info("MAFFT fallback ran on %d clusters in %.1fs", n_fallback, t_mafft)

    # Build full record set: include EVERY LTR member.
    refined_keys = set()
    all_records: List[RefRecord] = []
    for recs in cluster_records.values():
        for r in recs:
            all_records.append(r)
            refined_keys.add((r.chrom, r.start_orig, r.end_orig, r.role))
    for m in members:
        if (m.chrom, m.start, m.end, m.role) in refined_keys:
            continue
        all_records.append(RefRecord(
            chrom=m.chrom, start_orig=m.start, end_orig=m.end,
            strand=m.strand, role=m.role, parent=m.parent,
            lineage=m.lineage,
            cluster_id="", cluster_size=0,
            refinement_side="5" if m.role == "5LTR" else "3",
        ))

    # Per-element TSD recheck + revert.  Always populates tsd_outcome /
    # orig_tsd_*; reverts only if --no_tsd_revert was NOT passed.
    t0 = time.time()
    apply_tsd_gate_per_element(all_records, te_meta, genome, genome_lens,
                                enable_revert=not args.no_tsd_revert)
    t_tsd = time.time() - t0
    logger.info("TSD recheck + gate applied in %.2fs", t_tsd)

    # Outputs
    refined_gff = out_prefix + "_refined.gff3"
    per_element = out_prefix + "_per_element.tsv"
    cluster_tsv = out_prefix + "_clusters.tsv"
    run_json = out_prefix + "_run.json"

    counts = emit_refined_gff3(raw_lines, all_records, refined_gff)
    emit_per_element_tsv(all_records, per_element)
    emit_cluster_manifest(all_records, cluster_index, cluster_tsv)

    # Aggregate stats for run.json
    n_total = len(all_records)
    n_inner = sum(1 for r in all_records if r.refinement_method == "parasail_inner")
    n_outer = sum(1 for r in all_records if r.refinement_method == "parasail_outer")
    n_mafft = sum(1 for r in all_records if r.refinement_method == "mafft")
    n_unref = sum(1 for r in all_records if r.refinement_method == "none")
    n_dual = sum(1 for r in all_records if r.confidence == "dual")
    n_diverg = sum(1 for r in all_records if r.confidence == "divergent")
    n_io = sum(1 for r in all_records if r.confidence == "inner_only")
    n_reverted = sum(1 for r in all_records if r.revert_reason)
    tsd_counts = {k: 0 for k in
                  ("kept_exact", "kept_fuzzy", "shifted", "lost",
                   "gained", "both_none")}
    for r in all_records:
        if r.tsd_outcome and r.tsd_outcome in tsd_counts:
            tsd_counts[r.tsd_outcome] += 1

    summary = {
        "version": __SCHEMA_VERSION__,
        "policy": "inner_primary",
        "params": {
            "gff3": os.path.abspath(args.gff3),
            "genome": os.path.abspath(args.genome),
            "output": os.path.abspath(out_prefix),
            "identity": args.identity,
            "min_cluster_size": args.min_cluster_size,
            "anchor_len": args.anchor_len,
            "flank_len": args.flank_len,
            "snap_window": args.snap_window,
            "no_tsd_revert": bool(args.no_tsd_revert),
            "mafft_fallback": args.mafft_fallback,
            "mafft_fallback_threshold": args.mafft_fallback_threshold,
            "threads": args.threads,
            "workers": args.workers,
        },
        "counts": {
            "n_ltr_features": n_total,
            "n_refined_parasail_inner": n_inner,
            "n_refined_parasail_outer": n_outer,
            "n_refined_mafft": n_mafft,
            "n_unrefined": n_unref,
            "n_reverted": n_reverted,
            "n_confidence_dual": n_dual,
            "n_confidence_divergent": n_diverg,
            "n_confidence_inner_only": n_io,
            "ltr_refined_in_gff": counts["ltr_refined"],
            "te_refined_in_gff": counts["te_refined"],
            "tsd_emitted_rows": counts["tsd_emitted"],
            "tsd_dropped_rows": counts["tsd_dropped"],
            "tsd_outcome": tsd_counts,
            "n_clusters_total": len(cluster_members),
            "n_clusters_qualifying": len(qualifying),
            "n_clusters_mafft_fallback": n_fallback,
        },
        "timing_s": {
            "cluster": t_cluster,
            "parasail": t_parasail,
            "mafft_fallback": t_mafft,
            "tsd_gate": t_tsd,
            "total": time.time() - t0_total,
        },
    }
    with open(run_json, "w") as f:
        json.dump(summary, f, indent=2)

    logger.info("Wrote refined GFF3:    %s", refined_gff)
    logger.info("Wrote per-element TSV: %s", per_element)
    logger.info("Wrote cluster TSV:     %s", cluster_tsv)
    logger.info("Wrote run summary:     %s", run_json)
    logger.info("Refinement summary: %d/%d LTR features refined "
                "(parasail_inner=%d, mafft=%d), %d reverted; "
                "confidence dual=%d divergent=%d inner_only=%d",
                n_inner + n_mafft, n_total, n_inner, n_mafft, n_reverted,
                n_dual, n_diverg, n_io)
    return summary


__SCHEMA_VERSION__ = 2


# ============================================================
# CLI
# ============================================================

def build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="dante_ltr_refine",
        description="Per-element LTR boundary refinement (v2 — inner-primary "
                    "policy with TSD-aware gating).  See "
                    "docs/refine_v2_analysis.md for the design.",
    )
    ap.add_argument("-g", "--gff3", required=True,
                    help="DANTE_LTR GFF3 input")
    ap.add_argument("-s", "--genome", required=True,
                    help="Genome FASTA")
    ap.add_argument("-o", "--output", required=True,
                    help="Output prefix (files <prefix>_refined.gff3, "
                         "<prefix>_per_element.tsv, <prefix>_clusters.tsv, "
                         "<prefix>_run.json will be written)")
    ap.add_argument("--identity", type=float, default=0.9,
                    help="MMseqs2 cluster identity (default: %(default)s)")
    ap.add_argument("--min_cluster_size", type=int, default=6,
                    help="Min cluster members per role; cluster qualifies "
                         "when N(5'LTR) >= N AND N(3'LTR) >= N "
                         "(default: %(default)s)")
    ap.add_argument("--anchor_len", type=int, default=50,
                    help="Anchor length inside LTR, bp (default: %(default)s)")
    ap.add_argument("--flank_len", type=int, default=1000,
                    help="Flank length for parasail extension, bp "
                         "(default: %(default)s)")
    ap.add_argument("--snap_window", type=int, default=5,
                    help="±bp window for TG/CA motif snap (default: %(default)s)")
    ap.add_argument("--no_tsd_revert", action="store_true",
                    help="Disable the per-element TSD-loss revert rule "
                         "(diagnostic).  Refined coords are emitted even "
                         "when they destroy the original TSD.")
    ap.add_argument("--mafft_fallback_threshold", type=float, default=0.5,
                    help="Cluster-level inner-pool validation rate below "
                         "which MAFFT fallback fires (default: %(default)s)")
    ap.add_argument("--no-mafft-fallback", dest="mafft_fallback",
                    action="store_false",
                    help="Disable MAFFT fallback (parasail-only)")
    ap.set_defaults(mafft_fallback=True)
    ap.add_argument("--mafft_flank", type=int, default=50,
                    help="±bp flank for MAFFT fallback (default: %(default)s)")
    ap.add_argument("--mafft_timeout", type=int, default=600,
                    help="MAFFT fallback per-cluster timeout, s "
                         "(default: %(default)s)")

    # Parasail scoring (rarely changed)
    ap.add_argument("--match", type=int, default=2)
    ap.add_argument("--mismatch", type=int, default=-2)
    ap.add_argument("--gap_open", type=int, default=12)
    ap.add_argument("--gap_extend", type=int, default=3)
    ap.add_argument("--score_threshold", type=int, default=20)
    ap.add_argument("--min_num_alignments", type=int, default=3,
                    help="Minimum number of pairwise alignments per member "
                         "(Nth-largest aggregation; default: %(default)s)")
    ap.add_argument("--max_pairs", type=int, default=2000,
                    help="Cap pairs per cluster (default: %(default)s; 0=no cap)")
    ap.add_argument("--threads", type=int, default=4,
                    help="Threads for MMseqs2 / MAFFT (default: %(default)s)")
    ap.add_argument("--workers", type=int, default=4,
                    help="Parallel cluster workers (default: %(default)s)")
    ap.add_argument("-v", "--verbose", action="store_true")
    return ap


def configure_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
    logger.handlers = [handler]
    logger.setLevel(level)
    logger.propagate = False


def main(argv: Optional[List[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    configure_logging(args.verbose)
    refine_all(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
