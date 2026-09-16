#!/usr/bin/env python3
"""Measure core-mode constraints from existing dante_ltr annotations.

Core-domain mode (docs/core_domain_mode_design.md) seeds on the ordered
RT/RH/INT core and then searches outward for the LTR.  How far outward,
and how far apart the core domains may sit, are empirical quantities.
This script measures them on genomes that lineage mode has already
annotated, where both the LTR boundaries and the domain positions are
known, and writes databases/core_domain_order.csv.

Only elements of rank DLT/DLTP are used: a target site duplication was
found, so the boundaries are trustworthy.

Input is one or more dante_ltr output GFF3 files.  Standard library
only; the files are streamed, so whole-genome annotations are cheap.

Usage:
    utils/calibrate_core_constraints.py -g 'runs/*/dante_ltr*.gff3' \\
        -o databases/core_domain_order.csv \\
        --report core_constraints_report.tsv
"""

import argparse
import glob
import math
import os
import re
import sys
from collections import Counter, defaultdict

CORE = ('RT', 'RH', 'INT')

# Core order in element orientation, per superfamily.  This is the
# premise core mode rests on; the report measures how often it holds.
CORE_ORDER = {
    'Class_I/LTR/Ty1_copia': ('INT', 'RT', 'RH'),
    'Class_I/LTR/Ty3_gypsy': ('RT', 'RH', 'INT'),
    }

# Rank of an element whose boundaries we trust for measurement.
TRUSTED_RANKS = ('DLT', 'DLTP')


def superfamily_of(classification):
    """Map a REXdb classification onto one of the two LTR superfamilies.

    Returns None for anything that is not an LTR retrotransposon.
    """
    if not classification:
        return None
    parts = classification.replace('|', '/').split('/')
    # expect Class_I / LTR / Ty1_copia|Ty3_gypsy / ...
    if len(parts) < 3 or parts[0] != 'Class_I' or parts[1] != 'LTR':
        return None
    sf = parts[2].replace('/', '_')
    if sf.startswith('Ty1'):
        return 'Class_I/LTR/Ty1_copia'
    if sf.startswith('Ty3'):
        return 'Class_I/LTR/Ty3_gypsy'
    return None


def parse_attributes(field):
    """GFF3 column 9 -> dict.  Values are left URL-encoded as-is."""
    out = {}
    for item in field.rstrip(';').split(';'):
        if not item:
            continue
        key, _, value = item.partition('=')
        out[key.strip()] = value
    return out


def iter_gff3(path):
    """Yield (seqid, type, start, end, strand, attributes) per feature."""
    with open(path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            col = line.rstrip('\n').split('\t')
            if len(col) < 9:
                continue
            yield (col[0], col[2], int(col[3]), int(col[4]), col[6],
                   parse_attributes(col[8]))


class Element(object):
    """One transposable_element and the features that point at it."""

    __slots__ = ('eid', 'start', 'end', 'strand', 'rank', 'classification',
                 'ltr5', 'ltr3', 'domains')

    def __init__(self, eid, start, end, strand, rank, classification):
        self.eid = eid
        self.start = start
        self.end = end
        self.strand = strand
        self.rank = rank
        self.classification = classification
        self.ltr5 = None          # (start, end)
        self.ltr3 = None
        self.domains = []         # (start, end, name, classification)


def read_elements(path):
    """Collect elements and their children from one GFF3 file.

    Two passes over the file rather than one pass plus a buffer of every
    child feature: whole-genome annotations reach the gigabyte range and
    only the elements worth measuring are ever held in memory.
    """
    elements = {}
    for seqid, ftype, start, end, strand, attr in iter_gff3(path):
        if ftype != 'transposable_element':
            continue
        eid = attr.get('ID')
        # only ranks we trust are ever measured -- drop the rest now so
        # their children are skipped in the second pass too
        if eid is None or attr.get('Rank', '') not in TRUSTED_RANKS:
            continue
        elements[eid] = Element(eid, start, end, strand,
                                attr.get('Rank', ''),
                                attr.get('Final_Classification', ''))

    for seqid, ftype, start, end, strand, attr in iter_gff3(path):
        if ftype not in ('long_terminal_repeat', 'protein_domain'):
            continue
        el = elements.get(attr.get('Parent'))
        if el is None:
            continue
        if ftype == 'long_terminal_repeat':
            if attr.get('LTR') == '5LTR':
                el.ltr5 = (start, end)
            elif attr.get('LTR') == '3LTR':
                el.ltr3 = (start, end)
        else:
            el.domains.append((start, end, attr.get('Name', ''),
                               attr.get('Final_Classification', ''),
                               attr.get('Best_Hit_DB_Pos', ''),
                               attr.get('Similarity', ''),
                               attr.get('Relat_Length', '')))
    return list(elements.values())


DB_POS_RE = re.compile(r'^(\d+):(\d+)of(\d+)$')


def parse_db_pos(value):
    """'1:120of121' -> (1, 120, 121).  None when absent or malformed."""
    m = DB_POS_RE.match(value or '')
    if not m:
        return None
    return int(m.group(1)), int(m.group(2)), int(m.group(3))


def as_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def classify_pair(a, b, plus, args):
    """Is (a, b) one domain split in two, or two separate domains?

    A frameshift or a stop codon can break a single protein domain into
    two DANTE annotations.  The signature is that the two hits tile
    *complementary* parts of the same reference domain and sit next to
    each other on the genome.  Two genuine copies instead hit overlapping
    parts of the reference.

    Returns a record; 'split' is the verdict.
    """
    gap = (b[0] - a[1] - 1) if plus else (a[0] - b[1] - 1)
    pa, pb = parse_db_pos(a[4]), parse_db_pos(b[4])
    rec = {
        'name': a[2],
        'gap': gap,
        'db_overlap_frac': None,
        'db_ordered': None,
        'relat_length_sum': None,
        'split': False,
        }
    ra, rb = as_float(a[6]), as_float(b[6])
    if ra is not None and rb is not None:
        rec['relat_length_sum'] = ra + rb
    if pa is None or pb is None:
        return rec
    # overlap of the two reference-domain intervals, relative to the
    # shorter of the two
    overlap = min(pa[1], pb[1]) - max(pa[0], pb[0]) + 1
    shorter = min(pa[1] - pa[0] + 1, pb[1] - pb[0] + 1)
    rec['db_overlap_frac'] = max(0, overlap) / float(shorter) if shorter else 1.0
    # fragments of one domain appear in reference order along the element
    rec['db_ordered'] = pb[0] >= pa[0]
    rec['split'] = (rec['db_ordered'] and
                    rec['db_overlap_frac'] <= args.split_max_db_overlap and
                    gap <= args.split_max_gap)
    return rec


def merge_fragments(domains, plus, args):
    """Collapse split annotations of one domain into single entries.

    Returns (merged, pairs) where merged entries are
    (start, end, name, classification, n_fragments) and pairs holds a
    record for every same-name adjacent pair examined.
    """
    merged = []
    pairs = []
    for dom in domains:
        if merged and merged[-1][2] == dom[2]:
            prev_raw = merged[-1][5]
            pair = classify_pair(prev_raw, dom, plus, args)
            pairs.append(pair)
            if pair['split']:
                lo = min(merged[-1][0], dom[0])
                hi = max(merged[-1][1], dom[1])
                merged[-1] = (lo, hi, dom[2], merged[-1][3],
                              merged[-1][4] + 1, dom)
                continue
        merged.append((dom[0], dom[1], dom[2], dom[3], 1, dom))
    return merged, pairs


def measure(el, args):
    """Measure one element.  Returns a dict, or None if unusable.

    All distances are reported in element orientation, so a minus-strand
    element reads the same way as a plus-strand one.
    """
    if el.ltr5 is None or el.ltr3 is None or not el.domains:
        return None
    sf = superfamily_of(el.classification)
    if sf is None:
        return None

    plus = el.strand == '+'
    # sort domains in element orientation (5' first)
    domains = sorted(el.domains, key=lambda d: d[0], reverse=not plus)

    raw_core_names = [d[2] for d in domains if d[2] in CORE]
    raw_counts = Counter(raw_core_names)

    merged, pairs = merge_fragments(domains, plus, args)
    core = [d for d in merged if d[2] in CORE]
    core_names = tuple(d[2] for d in core)
    counts = Counter(core_names)

    rec = {
        'superfamily': sf,
        'rank': el.rank,
        'te_length': el.end - el.start + 1,
        'ltr5_width': el.ltr5[1] - el.ltr5[0] + 1,
        'ltr3_width': el.ltr3[1] - el.ltr3[0] + 1,
        # three levels of "has a core", see the report
        'core_types_present': set(raw_counts) == set(CORE),
        'core_clean_raw': (set(raw_counts) == set(CORE) and
                           all(raw_counts[n] == 1 for n in CORE)),
        'core_clean_merged': (set(counts) == set(CORE) and
                              all(counts[n] == 1 for n in CORE)),
        'core_multi_raw': any(raw_counts[n] > 1 for n in CORE),
        'core_recovered_by_merge': False,
        'pairs': pairs,
        'similarity': [(d[2], as_float(d[5]), as_float(d[6]))
                       for d in domains],
        'order_ok': None,
        'gaps': [],
        'core_span': None,
        'offset5prime': None,
        'offset3prime': None,
        'accessory_5': [],
        'accessory_3': [],
        'accessory_sf_mismatch': 0,
        'accessory_repeat_raw': False,
        }

    rec['core_recovered_by_merge'] = (rec['core_clean_merged'] and
                                      not rec['core_clean_raw'])

    # a repeated accessory domain would trip transparency rule 4 of the
    # design's blocking walk -- count it before merging, and after
    acc_counts_raw = Counter(d[2] for d in domains if d[2] not in CORE)
    acc_counts = Counter(d[2] for d in merged if d[2] not in CORE)
    rec['accessory_repeat_raw'] = any(c > 1 for c in acc_counts_raw.values())
    rec['accessory_repeat_merged'] = any(c > 1 for c in acc_counts.values())

    if not rec['core_clean_merged']:
        return rec

    rec['order_ok'] = core_names == CORE_ORDER[sf]

    # gaps between consecutive core domains, element orientation
    for prev, nxt in zip(core, core[1:]):
        gap = (nxt[0] - prev[1]) if plus else (prev[0] - nxt[1])
        rec['gaps'].append(gap - 1)

    core_lo = min(d[0] for d in core)
    core_hi = max(d[1] for d in core)
    rec['core_span'] = core_hi - core_lo + 1

    # distance from the element's 5' terminus to the 5'-most core base,
    # and from the 3'-most core base to the element's 3' terminus
    if plus:
        rec['offset5prime'] = core_lo - el.ltr5[0]
        rec['offset3prime'] = el.ltr3[1] - core_hi
    else:
        rec['offset5prime'] = el.ltr5[1] - core_hi
        rec['offset3prime'] = core_lo - el.ltr3[0]

    # accessory domains, by side of the core in element orientation
    for entry in merged:
        start, end, name, cls = entry[0], entry[1], entry[2], entry[3]
        if name in CORE:
            continue
        side5 = (end < core_lo) if plus else (start > core_hi)
        (rec['accessory_5'] if side5 else rec['accessory_3']).append(name)
        dsf = superfamily_of(cls)
        if dsf is not None and dsf != sf:
            rec['accessory_sf_mismatch'] += 1

    return rec


def quantile(values, q):
    """Linear-interpolation quantile of a list of numbers."""
    if not values:
        return None
    ordered = sorted(values)
    if len(ordered) == 1:
        return float(ordered[0])
    pos = (len(ordered) - 1) * q
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return float(ordered[lo])
    return ordered[lo] + (ordered[hi] - ordered[lo]) * (pos - lo)


def round_up(value, step):
    return int(math.ceil(value / float(step)) * step)


def round_down(value, step):
    return int(math.floor(value / float(step)) * step)


def _similarity_quantiles(recs):
    """Similarity and Relat_Length spread, core domains vs the rest.

    Both are measured on domains that already passed dante_filtering, so
    the low tail shows how close real core domains come to the current
    threshold -- i.e. how much a lower core-specific threshold could buy.
    """
    out = {}
    for group, wanted in (('core', True), ('accessory', False)):
        sims, rels = [], []
        for rec in recs:
            for name, sim, rel in rec['similarity']:
                if (name in CORE) != wanted:
                    continue
                if sim is not None:
                    sims.append(sim)
                if rel is not None:
                    rels.append(rel)
        out[group] = {
            'similarity': [quantile(sims, q)
                           for q in (0.0, 0.01, 0.05, 0.5)],
            'relat_length': [quantile(rels, q)
                             for q in (0.0, 0.01, 0.05, 0.5)],
            }
    return out


def summarize(records, args):
    """Group measurements by superfamily and derive the table rows."""
    by_sf = defaultdict(list)
    for rec in records:
        by_sf[rec['superfamily']].append(rec)

    rows = {}
    diagnostics = {}
    for sf, recs in sorted(by_sf.items()):
        usable = [r for r in recs if r['core_clean_merged']]
        geom = [r for r in usable if r['core_span'] is not None]

        gaps = [g for r in geom for g in r['gaps']]
        spans = [r['core_span'] for r in geom]
        off5 = [r['offset5prime'] for r in geom]
        off3 = [r['offset3prime'] for r in geom]
        ltrs = [w for r in geom for w in (r['ltr5_width'], r['ltr3_width'])]
        lengths = [r['te_length'] for r in geom]

        acc5 = defaultdict(int)
        acc3 = defaultdict(int)
        for r in geom:
            for name in set(r['accessory_5']):
                acc5[name] += 1
            for name in set(r['accessory_3']):
                acc3[name] += 1

        n_geom = max(len(geom), 1)
        keep5 = sorted((n for n, c in acc5.items()
                        if c / n_geom >= args.accessory_min_freq),
                       key=lambda n: -acc5[n])
        keep3 = sorted((n for n, c in acc3.items()
                        if c / n_geom >= args.accessory_min_freq),
                       key=lambda n: -acc3[n])

        hi = args.upper_quantile
        m = args.margin
        rows[sf] = {
            'Superfamily': sf,
            'core_order': ' '.join(CORE_ORDER[sf]),
            'offset5prime': round_up(m * (quantile(off5, hi) or 0), 500),
            'offset3prime': round_up(m * (quantile(off3, hi) or 0), 500),
            'core_max_gap': round_up(m * (quantile(gaps, hi) or 0), 500),
            'core_max_span': round_up(m * (quantile(spans, hi) or 0), 500),
            'ltr_length': max(args.min_ltr_floor,
                              round_down(quantile(ltrs,
                                                  args.lower_quantile) or 0,
                                         10)),
            'accessory_5': ' '.join(keep5) if keep5 else '-',
            'accessory_3': ' '.join(keep3) if keep3 else '-',
            }

        order_known = [r for r in usable if r['order_ok'] is not None]
        n = float(len(recs)) if recs else 1.0
        split_pairs = [p for r in recs for p in r['pairs'] if p['split']]
        dup_pairs = [p for r in recs for p in r['pairs'] if not p['split']]
        diagnostics[sf] = {
            'elements_trusted_rank': len(recs),
            'core_types_present': sum(1 for r in recs
                                      if r['core_types_present']),
            'core_clean_raw': sum(1 for r in recs if r['core_clean_raw']),
            'core_complete': len(usable),
            'core_complete_frac': len(usable) / n,
            'core_clean_raw_frac': sum(1 for r in recs
                                       if r['core_clean_raw']) / n,
            'core_recovered_by_merge': sum(1 for r in recs
                                           if r['core_recovered_by_merge']),
            'core_multi_raw': sum(1 for r in recs if r['core_multi_raw']),
            'accessory_repeat_raw': sum(1 for r in recs
                                        if r['accessory_repeat_raw']),
            'accessory_repeat_merged': sum(1 for r in recs
                                           if r['accessory_repeat_merged']),
            'split_pairs': len(split_pairs),
            'dup_pairs': len(dup_pairs),
            'split_by_name': Counter(p['name'] for p in split_pairs),
            'dup_by_name': Counter(p['name'] for p in dup_pairs),
            'split_gap_q': [quantile([p['gap'] for p in split_pairs], q)
                            for q in (0.5, 0.95, 0.99, 1.0)],
            'split_relat_q': [quantile([p['relat_length_sum']
                                        for p in split_pairs
                                        if p['relat_length_sum'] is not None],
                                       q)
                              for q in (0.05, 0.5, 0.95, 1.0)],
            'similarity_q': _similarity_quantiles(recs),
            'order_concordant': sum(1 for r in order_known if r['order_ok']),
            'order_discordant': sum(1 for r in order_known
                                    if not r['order_ok']),
            'accessory_sf_mismatch': sum(r['accessory_sf_mismatch']
                                         for r in geom),
            'quantiles': {
                'gap': [quantile(gaps, q) for q in (0.5, 0.95, 0.99, 1.0)],
                'core_span': [quantile(spans, q)
                              for q in (0.5, 0.95, 0.99, 1.0)],
                'offset5prime': [quantile(off5, q)
                                 for q in (0.5, 0.95, 0.99, 1.0)],
                'offset3prime': [quantile(off3, q)
                                 for q in (0.5, 0.95, 0.99, 1.0)],
                'ltr_width': [quantile(ltrs, q)
                              for q in (0.01, 0.5, 0.95, 1.0)],
                'te_length': [quantile(lengths, q)
                              for q in (0.5, 0.95, 0.99, 1.0)],
                },
            }
    return rows, diagnostics


COLUMNS = ('Superfamily', 'core_order', 'offset5prime', 'offset3prime',
           'core_max_gap', 'core_max_span', 'ltr_length', 'accessory_5',
           'accessory_3')


def write_table(rows, path):
    with open(path, 'w') as fh:
        fh.write('\t'.join(COLUMNS) + '\n')
        for sf in sorted(rows):
            fh.write('\t'.join(str(rows[sf][c]) for c in COLUMNS) + '\n')


def write_report(diagnostics, path):
    labels = ('q50', 'q95', 'q99', 'max')
    with open(path, 'w') as fh:
        fh.write('superfamily\tmeasure\t' + '\t'.join(labels) + '\n')
        for sf in sorted(diagnostics):
            for measure, values in sorted(
                    diagnostics[sf]['quantiles'].items()):
                cells = ['' if v is None else '{:.0f}'.format(v)
                         for v in values]
                fh.write('{}\t{}\t{}\n'.format(sf, measure,
                                               '\t'.join(cells)))


def print_diagnostics(diagnostics, out=sys.stdout):
    for sf in sorted(diagnostics):
        d = diagnostics[sf]
        out.write('\n{}\n'.format(sf))
        out.write('  elements at rank DLT/DLTP        : {}\n'
                  .format(d['elements_trusted_rank']))
        out.write('  all three core types present     : {}\n'
                  .format(d['core_types_present']))
        out.write('  ...exactly one annotation each   : {} ({:.1%})\n'
                  .format(d['core_clean_raw'], d['core_clean_raw_frac']))
        out.write('  ...after merging split domains   : {} ({:.1%})'
                  '  [+{} recovered]\n'
                  .format(d['core_complete'], d['core_complete_frac'],
                          d['core_recovered_by_merge']))
        out.write('  elements with a repeated core domain    : {}\n'
                  .format(d['core_multi_raw']))
        out.write('  elements with a repeated accessory dom. : {} '
                  '-> {} after merging\n'
                  .format(d['accessory_repeat_raw'],
                          d['accessory_repeat_merged']))
        out.write('  same-name adjacent pairs: {} split / {} duplicate\n'
                  .format(d['split_pairs'], d['dup_pairs']))
        if d['split_by_name']:
            out.write('    split by domain    : {}\n'.format(
                ', '.join('{} {}'.format(k, v) for k, v in
                          d['split_by_name'].most_common())))
        if d['dup_by_name']:
            out.write('    duplicate by domain: {}\n'.format(
                ', '.join('{} {}'.format(k, v) for k, v in
                          d['dup_by_name'].most_common())))
        fmt = lambda vals: ''.join(
            '{:>9}'.format('-' if v is None else '{:.2f}'.format(v))
            for v in vals)
        out.write('    split gap bp (q50 q95 q99 max) : {}\n'.format(
            ''.join('{:>8}'.format('-' if v is None else '{:.0f}'.format(v))
                    for v in d['split_gap_q'])))
        out.write('    split summed Relat_Length      : {}   '
                  '(q05 q50 q95 max)\n'.format(fmt(d['split_relat_q'])))
        sq = d['similarity_q']
        out.write('  domain quality, min/q01/q05/q50 (post-filter)\n')
        for group in ('core', 'accessory'):
            out.write('    {:<10} Similarity   : {}\n'
                      .format(group, fmt(sq[group]['similarity'])))
            out.write('    {:<10} Relat_Length : {}\n'
                      .format(group, fmt(sq[group]['relat_length'])))
        total_order = d['order_concordant'] + d['order_discordant']
        if total_order:
            out.write('  core order matches superfamily   : {} / {} '
                      '({} discordant)\n'
                      .format(d['order_concordant'], total_order,
                              d['order_discordant']))
        out.write('  accessory domains whose superfamily disagrees: {}\n'
                  .format(d['accessory_sf_mismatch']))
        q = d['quantiles']
        out.write('  {:<14}{:>9}{:>9}{:>9}{:>9}\n'
                  .format('', 'q50', 'q95', 'q99', 'max'))
        for measure in ('gap', 'core_span', 'offset5prime', 'offset3prime',
                        'te_length'):
            vals = ''.join('{:>9}'.format('-' if v is None
                                          else '{:.0f}'.format(v))
                           for v in q[measure])
            out.write('  {:<14}{}\n'.format(measure, vals))
        vals = ''.join('{:>9}'.format('-' if v is None
                                      else '{:.0f}'.format(v))
                       for v in q['ltr_width'])
        out.write('  {:<14}{}   (q01 q50 q95 max)\n'.format('ltr_width',
                                                            vals))


def get_arguments():
    parser = argparse.ArgumentParser(
        description=('Measure core-mode constraints (core gap/span, LTR '
                     'search offsets, accessory domain inventory) from '
                     'existing dante_ltr annotations, and write '
                     'core_domain_order.csv.'),
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        '-g', '--gff3', required=True, nargs='+', action='store',
        help=('dante_ltr output GFF3 file(s).  Shell globs are accepted '
              'and are also expanded internally, so quoted patterns such '
              "as 'runs/*/dante_ltr*.gff3' work."))
    parser.add_argument(
        '-o', '--output', default=None,
        help='write the derived constraints table here (TSV)')
    parser.add_argument(
        '--report', default=None,
        help='write the full quantile report here (TSV)')
    parser.add_argument(
        '--upper_quantile', type=float, default=0.99,
        help=('quantile used for the cap-like values -- offsets, gap and '
              'span (default %(default)s)'))
    parser.add_argument(
        '--lower_quantile', type=float, default=0.01,
        help=('quantile used for ltr_length, which is a minimum-length '
              'filter rather than a cap (default %(default)s)'))
    parser.add_argument(
        '--margin', type=float, default=1.0,
        help=('multiply the cap-like values by this before rounding.  The\n'
              'measurement is censored -- it can only see elements that\n'
              'lineage mode already found, i.e. elements whose LTR fell\n'
              'inside lineage mode\'s own search window -- so the observed\n'
              'offsets are a lower bound on what core mode needs.  Use\n'
              '1.5 when deriving the shipped table (default %(default)s)'))
    parser.add_argument(
        '--split_max_gap', type=int, default=1000,
        help=('two same-name annotations closer than this may be one\n'
              'domain split by a frameshift (default %(default)s)'))
    parser.add_argument(
        '--split_max_db_overlap', type=float, default=0.3,
        help=('...and only if their hits to the reference domain overlap\n'
              'by at most this fraction of the shorter hit -- genuine\n'
              'duplicates hit the same part twice (default %(default)s)'))
    parser.add_argument(
        '--min_ltr_floor', type=int, default=100,
        help='ltr_length is never set below this (default %(default)s)')
    parser.add_argument(
        '--accessory_min_freq', type=float, default=0.02,
        help=('an accessory domain is listed for a side when it occurs '
              'there in at least this fraction of elements '
              '(default %(default)s)'))
    return parser.parse_args()


def main():
    args = get_arguments()

    paths = []
    for pattern in args.gff3:
        matched = sorted(glob.glob(pattern))
        paths.extend(matched if matched else [pattern])
    missing = [p for p in paths if not os.path.isfile(p)]
    if missing:
        sys.stderr.write('ERROR: not found: {}\n'.format(', '.join(missing)))
        sys.exit(1)
    if not paths:
        sys.stderr.write('ERROR: no input files matched\n')
        sys.exit(1)

    records = []
    for path in paths:
        elements = read_elements(path)
        measured = [m for m in (measure(el, args) for el in elements)
                    if m is not None]
        records.extend(measured)
        print('{}: {} elements, {} usable at rank DLT/DLTP'
              .format(path, len(elements), len(measured)))

    if not records:
        sys.stderr.write('ERROR: no rank DLT/DLTP LTR elements found; '
                         'nothing to measure\n')
        sys.exit(1)

    rows, diagnostics = summarize(records, args)
    print_diagnostics(diagnostics)

    if args.output:
        write_table(rows, args.output)
        print('\nconstraints table -> {}'.format(args.output))
    else:
        print('')
        print('\t'.join(COLUMNS))
        for sf in sorted(rows):
            print('\t'.join(str(rows[sf][c]) for c in COLUMNS))

    if args.report:
        write_report(diagnostics, args.report)
        print('quantile report   -> {}'.format(args.report))


if __name__ == '__main__':
    main()
