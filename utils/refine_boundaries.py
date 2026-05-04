"""Per-element LTR boundary refinement engine.

Reads a DANTE_LTR GFF3 + reference FASTA, runs MMseqs2 per-lineage
clustering, refines each cluster member's outward-facing LTR boundary
by parasail anchored extension (with optional MAFFT change-point
fallback for low-validation-rate clusters), and emits:

  * <prefix>_refined.gff3            -- Output A
  * <prefix>_per_element.tsv         -- per-element refinement record
  * <prefix>_clusters.tsv            -- cluster manifest with stats
  * <prefix>_run.json                -- parameters + timing

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
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Local lib
HERE = os.path.dirname(os.path.realpath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from parasail_boundary import (  # type: ignore
    Member,
    SUPPORTED_MOTIFS,
    process_cluster_side,
    revcomp,
)


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
      te_meta  : dict keyed by TE id -> {classification, rank, ...}
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
                    te_meta[te_id] = {
                        'classification': _attr(attrs, 'Final_Classification') or "",
                        'rank': _attr(attrs, 'Rank') or "",
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
    chrom: str
    start_orig: int
    end_orig: int
    strand: str
    role: str
    parent: str
    lineage: str
    cluster_id: str = ""
    cluster_size: int = 0
    refinement_side: str = ""  # "5" for 5LTRs (refines bio-5'); "3" for 3LTRs

    parasail_n_pairs: int = 0
    parasail_corrected_g: Optional[int] = None
    parasail_motif_at: Optional[bool] = None
    parasail_snap_offset: int = 0
    parasail_motif_ok: Optional[bool] = None

    mafft_corrected_g: Optional[int] = None
    mafft_motif_ok: Optional[bool] = None

    final_corrected_g: Optional[int] = None
    final_method: str = "none"        # parasail | mafft | none
    final_confidence: str = "unrefined"  # high | medium | low | unrefined
    final_motif_ok: Optional[bool] = None


def _classify_confidence(method: str, motif_ok: Optional[bool],
                         snap_offset: int) -> str:
    """Refinement_Confidence policy.

    high   : parasail or MAFFT validated, no snap (offset 0)
    medium : validated after snap (small adjustment)
    low    : refinement returned a position but motif validation failed
             (or could not be evaluated)
    unrefined : no refinement available
    """
    if method == "none":
        return "unrefined"
    if motif_ok is True and snap_offset == 0:
        return "high"
    if motif_ok is True:
        return "medium"
    return "low"


# ============================================================
# Cluster orchestration
# ============================================================

def refine_one_cluster(cluster_members: List[Member], cluster_id: str,
                       cluster_lineage: str,
                       genome, genome_lens: Dict[str, int],
                       opts: argparse.Namespace
                       ) -> List[RefRecord]:
    """Run parasail on both sides of a single cluster. MAFFT fallback,
    if requested, is wired in `refine_all`.  Returns one RefRecord per
    LTR feature in the cluster.
    """
    # Parasail side=5 (5'LTRs only); side=3 (3'LTRs only)
    rows_5, _ = process_cluster_side(
        cluster_members, "5", genome, genome_lens,
        opts.anchor_len, opts.flank_len,
        match=opts.match, mismatch=opts.mismatch,
        gap_open=opts.gap_open, gap_extend=opts.gap_extend,
        score_threshold=opts.score_threshold,
        min_num_alignments=opts.min_num_alignments,
        snap_window=opts.snap_window,
        motif=opts.boundary_motif,
        max_pairs=opts.max_pairs,
    )
    rows_3, _ = process_cluster_side(
        cluster_members, "3", genome, genome_lens,
        opts.anchor_len, opts.flank_len,
        match=opts.match, mismatch=opts.mismatch,
        gap_open=opts.gap_open, gap_extend=opts.gap_extend,
        score_threshold=opts.score_threshold,
        min_num_alignments=opts.min_num_alignments,
        snap_window=opts.snap_window,
        motif=opts.boundary_motif,
        max_pairs=opts.max_pairs,
    )

    by_id: Dict[Tuple[str, int, int, str], Dict] = {}
    for r in rows_5 + rows_3:
        key = (r["chrom"], r["start"], r["end"], r["role"])
        by_id[key] = r

    records: List[RefRecord] = []
    for m in cluster_members:
        rec = RefRecord(
            chrom=m.chrom, start_orig=m.start, end_orig=m.end,
            strand=m.strand, role=m.role, parent=m.parent, lineage=m.lineage,
            cluster_id=cluster_id,
            cluster_size=len(cluster_members),
            refinement_side="5" if m.role == "5LTR" else "3",
        )
        r = by_id.get((m.chrom, m.start, m.end, m.role))
        if r is not None:
            rec.parasail_n_pairs = r["n_pairs"]
            rec.parasail_corrected_g = r["corrected_g"]
            rec.parasail_motif_at = r["motif_at_parasail"]
            rec.parasail_snap_offset = r["snap_offset"]
            rec.parasail_motif_ok = r["motif_ok"]
            if r["corrected_g"] is not None:
                rec.final_corrected_g = r["corrected_g"]
                rec.final_method = "parasail"
                rec.final_motif_ok = r["motif_ok"]
                rec.final_confidence = _classify_confidence(
                    "parasail", r["motif_ok"], r["snap_offset"])
        records.append(rec)
    return records


def cluster_validation_rate(records: List[RefRecord]) -> Tuple[float, int, int]:
    """Fraction of records (in this cluster) where parasail-validated.

    Members the parasail engine couldn't resolve count as un-validated
    in the denominator (matches plan §5.3 wording).
    Returns (rate, n_validated, n_total).
    """
    n_total = len(records)
    if n_total == 0:
        return 1.0, 0, 0
    n_validated = sum(1 for r in records if r.parasail_motif_ok is True)
    return n_validated / n_total, n_validated, n_total


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
               "--boundary_motif", opts.boundary_motif]
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

        # Merge into records.  Only override when parasail did NOT validate
        # AND MAFFT's call validated (or, if motif=none, MAFFT produced a
        # call when parasail did not).
        for rec in records:
            nm = (f"{rec.chrom}__{rec.start_orig}__{rec.end_orig}"
                  f"__{rec.role}__{rec.strand}")
            mafft = by_name.get(nm)
            if mafft is None:
                continue
            cg, mok = mafft
            rec.mafft_corrected_g = cg
            rec.mafft_motif_ok = mok
            parasail_validated = rec.parasail_motif_ok is True
            if parasail_validated:
                continue
            if cg is None:
                continue
            mafft_validated = mok is True
            if opts.boundary_motif == "none":
                # Both trivially "validate" under motif=none; only override
                # if parasail produced no position at all.
                if rec.parasail_corrected_g is not None:
                    continue
            elif not mafft_validated:
                # Parasail did not validate; mafft did not either -> keep
                # parasail (or unrefined if parasail had no position).
                continue
            # Use MAFFT
            rec.final_corrected_g = cg
            rec.final_method = "mafft"
            rec.final_motif_ok = mok
            rec.final_confidence = _classify_confidence(
                "mafft", mok, snap_offset=0)
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
    final_method != "none", else "no".
    """
    s, e = rec.start_orig, rec.end_orig
    if rec.final_method == "none" or rec.final_corrected_g is None:
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


def emit_refined_gff3(raw_lines: List[str], records: List[RefRecord],
                      out_path: str, motif: str) -> Dict[str, int]:
    """Emit refined GFF3 by replaying input lines and updating LTR/TE rows.

    Returns counts dict for logging.
    """
    # Index records by (chrom, start_orig, end_orig, role)
    rec_idx: Dict[Tuple[str, int, int, str], RefRecord] = {}
    for r in records:
        rec_idx[(r.chrom, r.start_orig, r.end_orig, r.role)] = r

    # Index by parent for TE-envelope updates
    rec_by_parent: Dict[str, Dict[str, RefRecord]] = defaultdict(dict)
    for r in records:
        rec_by_parent[r.parent][r.role] = r

    counts = {"te_total": 0, "te_refined": 0,
              "ltr_total": 0, "ltr_refined": 0}

    drop_keys = ("Original_Start", "Original_End",
                 "Refinement_Method", "Refinement_Confidence",
                 "Cluster_ID", "Cluster_Size", "TG_OK", "CA_OK")

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
                add_attrs = (
                    f"Original_Start={rec.start_orig};"
                    f"Original_End={rec.end_orig};"
                    f"Refinement_Method={rec.final_method};"
                    f"Refinement_Confidence={rec.final_confidence};"
                    f"Cluster_ID={rec.cluster_id};"
                    f"Cluster_Size={rec.cluster_size}"
                )
                if motif == "TG/CA":
                    if role == "5LTR":
                        tg = ("TRUE" if rec.final_motif_ok is True else
                              "FALSE" if rec.final_motif_ok is False else "NA")
                        add_attrs += f";TG_OK={tg};CA_OK=NA"
                    else:
                        ca = ("TRUE" if rec.final_motif_ok is True else
                              "FALSE" if rec.final_motif_ok is False else "NA")
                        add_attrs += f";TG_OK=NA;CA_OK={ca}"
                new_attrs = attrs_clean + ";" + add_attrs if attrs_clean else add_attrs
                out.write("\t".join([seqname, src, ftype, str(ns), str(ne),
                                     sc, strand, ph, new_attrs]) + "\n")
            elif ftype == 'transposable_element':
                counts["te_total"] += 1
                te_id = _attr(attrs, 'ID') or ''
                pair = rec_by_parent.get(te_id, {})
                rec5 = pair.get("5LTR"); rec3 = pair.get("3LTR")

                # Compute envelope from refined LTR coords, falling back
                # to original where a side wasn't refined.
                refined_any = False
                env_s, env_e = int(s1), int(e1)
                if rec5 is not None or rec3 is not None:
                    if rec5 is not None:
                        ns5, ne5, ap5 = _refined_ltr_coords(rec5)
                        if ap5 == "yes":
                            refined_any = True
                    else:
                        ns5 = int(s1); ne5 = int(e1)
                    if rec3 is not None:
                        ns3, ne3, ap3 = _refined_ltr_coords(rec3)
                        if ap3 == "yes":
                            refined_any = True
                    else:
                        ns3 = int(s1); ne3 = int(e1)
                    env_s = min(ns5, ns3, int(s1))
                    env_e = max(ne5, ne3, int(e1))
                    if rec5 is not None and rec3 is not None:
                        env_s = min(ns5, ns3)
                        env_e = max(ne5, ne3)

                if refined_any:
                    counts["te_refined"] += 1

                method_5 = rec5.final_method if rec5 is not None else "none"
                method_3 = rec3.final_method if rec3 is not None else "none"
                conf_5 = rec5.final_confidence if rec5 is not None else "unrefined"
                conf_3 = rec3.final_confidence if rec3 is not None else "unrefined"
                conf_rank = {"high": 3, "medium": 2, "low": 1, "unrefined": 0}
                te_conf = min((conf_5, conf_3), key=lambda c: conf_rank.get(c, 0))
                te_method = (method_5 if method_5 != "none" else
                             method_3 if method_3 != "none" else "none")
                cluster_id = (rec5.cluster_id if rec5 is not None else
                              rec3.cluster_id if rec3 is not None else "")
                cluster_size = (rec5.cluster_size if rec5 is not None else
                                rec3.cluster_size if rec3 is not None else 0)

                add_attrs = (
                    f"Original_Start={s1};"
                    f"Original_End={e1};"
                    f"Refinement_Method={te_method};"
                    f"Refinement_Confidence={te_conf};"
                    f"Cluster_ID={cluster_id};"
                    f"Cluster_Size={cluster_size}"
                )
                if motif == "TG/CA":
                    tg = "NA"
                    ca = "NA"
                    if rec5 is not None:
                        tg = ("TRUE" if rec5.final_motif_ok is True else
                              "FALSE" if rec5.final_motif_ok is False else "NA")
                    if rec3 is not None:
                        ca = ("TRUE" if rec3.final_motif_ok is True else
                              "FALSE" if rec3.final_motif_ok is False else "NA")
                    add_attrs += f";TG_OK={tg};CA_OK={ca}"
                new_attrs = attrs_clean + ";" + add_attrs if attrs_clean else add_attrs
                out.write("\t".join([seqname, src, ftype, str(env_s), str(env_e),
                                     sc, strand, ph, new_attrs]) + "\n")
            else:
                out.write(line)
    return counts


# ============================================================
# Per-element TSV
# ============================================================

def emit_per_element_tsv(records: List[RefRecord], out_path: str) -> None:
    cols = [
        "chrom", "start_orig", "end_orig", "strand", "role",
        "parent_te_id", "lineage_full",
        "cluster_id", "cluster_size", "refinement_side",
        "parasail_n_pairs", "parasail_corrected_g",
        "parasail_motif_ok", "parasail_snap_offset",
        "mafft_corrected_g", "mafft_motif_ok",
        "final_corrected_g", "final_method", "final_confidence",
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
            anchored = (r.start_orig if (r.role == "5LTR" and r.strand == '+') or
                                        (r.role == "3LTR" and r.strand == '-')
                        else r.end_orig)
            shift = (r.final_corrected_g - anchored
                     if r.final_corrected_g is not None else None)
            row = [
                r.chrom, r.start_orig, r.end_orig, r.strand, r.role,
                r.parent, r.lineage,
                r.cluster_id, r.cluster_size, r.refinement_side,
                r.parasail_n_pairs, r.parasail_corrected_g,
                r.parasail_motif_ok, r.parasail_snap_offset,
                r.mafft_corrected_g, r.mafft_motif_ok,
                r.final_corrected_g, r.final_method, r.final_confidence,
                shift,
            ]
            f.write("\t".join(fmt(v) for v in row) + "\n")


# ============================================================
# Cluster manifest
# ============================================================

def emit_cluster_manifest(records: List[RefRecord], cluster_index,
                          out_path: str) -> None:
    """cluster_index: dict cluster_id -> dict with 'lineage', 'mafft_invoked',
    'mafft_validation_rate'."""
    by_cl: Dict[str, List[RefRecord]] = defaultdict(list)
    for r in records:
        if r.cluster_id:
            by_cl[r.cluster_id].append(r)

    with open(out_path, "w") as f:
        f.write("cluster_id\tlineage_full\tn_total\tn_5ltr\tn_3ltr\t"
                "n_validated_5\tn_validated_3\t"
                "parasail_validation_rate\tmafft_invoked\t"
                "mafft_validation_rate\tlow_confidence\n")
        for cid, recs in by_cl.items():
            lin = recs[0].lineage if recs else ""
            n5 = sum(1 for r in recs if r.role == "5LTR")
            n3 = sum(1 for r in recs if r.role == "3LTR")
            v5 = sum(1 for r in recs if r.role == "5LTR" and r.parasail_motif_ok is True)
            v3 = sum(1 for r in recs if r.role == "3LTR" and r.parasail_motif_ok is True)
            denom = len(recs) or 1
            psrate = (sum(1 for r in recs if r.parasail_motif_ok is True) / denom)
            cmeta = cluster_index.get(cid, {})
            mafft_inv = "TRUE" if cmeta.get("mafft_invoked") else "FALSE"
            mafft_rate = cmeta.get("mafft_validation_rate")
            mafft_rate_str = f"{mafft_rate:.4f}" if mafft_rate is not None else "NA"
            low_conf = "TRUE" if (v5 + v3) < 4 else "FALSE"
            f.write(f"{cid}\t{lin}\t{len(recs)}\t{n5}\t{n3}\t"
                    f"{v5}\t{v3}\t{psrate:.4f}\t{mafft_inv}\t"
                    f"{mafft_rate_str}\t{low_conf}\n")


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

    # Run parasail per cluster (qualify by min_cluster_size)
    qualifying = [(ci, mems, cluster_lineage[ci], cluster_rep[ci])
                  for ci, mems in enumerate(cluster_members)
                  if len(mems) >= args.min_cluster_size]
    logger.info("  %d clusters qualify (>= --min_cluster_size=%d) for parasail",
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

    # MAFFT fallback per cluster
    t0 = time.time()
    n_fallback = 0
    if args.mafft_fallback and args.boundary_motif != "none":
        # Build cluster_id -> members map for MAFFT call
        cl_mems_by_id: Dict[str, List[Member]] = {}
        for ci, mems in enumerate(cluster_members):
            if len(mems) < args.min_cluster_size:
                continue
            cluster_id = f"clu_{ci:06d}_{cluster_rep[ci][:80]}"
            cl_mems_by_id[cluster_id] = mems
        for cid, recs in cluster_records.items():
            rate, _, _ = cluster_validation_rate(recs)
            if rate < args.mafft_fallback_threshold:
                n_fallback += 1
                cluster_index[cid]["mafft_invoked"] = True
                run_mafft_fallback(recs, cl_mems_by_id[cid],
                                   cid, cluster_index[cid]["lineage"],
                                   genome, genome_lens, args)
                # Re-rate after MAFFT merge
                n_validated = sum(1 for r in recs if (
                    r.final_motif_ok is True))
                cluster_index[cid]["mafft_validation_rate"] = (
                    n_validated / len(recs) if recs else 0.0)
    t_mafft = time.time() - t0
    logger.info("MAFFT fallback ran on %d clusters in %.1fs", n_fallback, t_mafft)

    # Build full record set: include EVERY LTR member; un-refined ones get
    # default RefRecord.  Records from qualifying clusters are already
    # populated; un-refined members get a stub.
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

    # Outputs
    refined_gff = out_prefix + "_refined.gff3"
    per_element = out_prefix + "_per_element.tsv"
    cluster_tsv = out_prefix + "_clusters.tsv"
    run_json = out_prefix + "_run.json"

    counts = emit_refined_gff3(raw_lines, all_records, refined_gff,
                               args.boundary_motif)
    emit_per_element_tsv(all_records, per_element)
    emit_cluster_manifest(all_records, cluster_index, cluster_tsv)

    # Aggregate stats for run.json
    n_total = len(all_records)
    n_parasail = sum(1 for r in all_records if r.final_method == "parasail")
    n_mafft = sum(1 for r in all_records if r.final_method == "mafft")
    n_unref = sum(1 for r in all_records if r.final_method == "none")
    n_validated = sum(1 for r in all_records if r.final_motif_ok is True)
    n_high = sum(1 for r in all_records if r.final_confidence == "high")
    n_med = sum(1 for r in all_records if r.final_confidence == "medium")
    n_low = sum(1 for r in all_records if r.final_confidence == "low")

    summary = {
        "version": __SCHEMA_VERSION__,
        "params": {
            "gff3": os.path.abspath(args.gff3),
            "genome": os.path.abspath(args.genome),
            "output": os.path.abspath(out_prefix),
            "identity": args.identity,
            "min_cluster_size": args.min_cluster_size,
            "anchor_len": args.anchor_len,
            "flank_len": args.flank_len,
            "snap_window": args.snap_window,
            "boundary_motif": args.boundary_motif,
            "mafft_fallback": args.mafft_fallback,
            "mafft_fallback_threshold": args.mafft_fallback_threshold,
            "threads": args.threads,
            "workers": args.workers,
        },
        "counts": {
            "n_ltr_features": n_total,
            "n_refined_parasail": n_parasail,
            "n_refined_mafft": n_mafft,
            "n_unrefined": n_unref,
            "n_validated": n_validated,
            "n_high": n_high, "n_medium": n_med, "n_low": n_low,
            "ltr_refined_in_gff": counts["ltr_refined"],
            "te_refined_in_gff": counts["te_refined"],
            "n_clusters_total": len(cluster_members),
            "n_clusters_qualifying": len(qualifying),
            "n_clusters_mafft_fallback": n_fallback,
        },
        "timing_s": {
            "cluster": t_cluster,
            "parasail": t_parasail,
            "mafft_fallback": t_mafft,
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
                "(parasail=%d, mafft=%d), %d validated",
                n_parasail + n_mafft, n_total, n_parasail, n_mafft, n_validated)
    return summary


__SCHEMA_VERSION__ = 1


# ============================================================
# CLI
# ============================================================

def build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="dante_ltr_refine",
        description="Per-element LTR boundary refinement (parasail anchored "
                    "extension + optional MAFFT change-point fallback).",
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
                    help="Min cluster members for parasail (default: %(default)s)")
    ap.add_argument("--anchor_len", type=int, default=50,
                    help="Anchor length inside LTR, bp (default: %(default)s)")
    ap.add_argument("--flank_len", type=int, default=1000,
                    help="Flank length for parasail extension, bp "
                         "(default: %(default)s)")
    ap.add_argument("--snap_window", type=int, default=5,
                    help="±bp window for motif snap (default: %(default)s)")
    ap.add_argument("--boundary_motif", default="TG/CA",
                    choices=list(SUPPORTED_MOTIFS),
                    help="Boundary motif policy (default: %(default)s)")
    ap.add_argument("--mafft_fallback_threshold", type=float, default=0.5,
                    help="Cluster-level parasail validation rate below which "
                         "MAFFT fallback fires (default: %(default)s)")
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
