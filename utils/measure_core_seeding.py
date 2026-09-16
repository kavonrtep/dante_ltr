#!/usr/bin/env python3
"""Measure core-mode seeding directly on raw DANTE output.

`calibrate_core_constraints.py` measures dante_ltr *output*, which can
only show elements lineage mode already accepted.  Two questions cannot
be answered there:

  * **Fragmented domains.**  A frameshift or stop codon can split one
    protein domain into two DANTE annotations.  Lineage mode discards
    any cluster with a repeated domain name outright
    (`clean_domain_clusters()`, ltr_utils.R:206 -- `N_unique_domains ==
    N_domains`), so affected elements never reach the output at all.
    They are invisible downstream, by construction.
  * **Filter thresholds.**  Domains that fail `dante_filtering` are
    absent from the output, so the effect of a lower threshold cannot be
    read off it.

Both need the raw DANTE GFF3.  This script applies a configurable
filter, optionally merges split annotations, enumerates ordered RT/RH/INT
core seeds per docs/core_domain_mode_design.md §6.2, and reports what
changes.

The GFF3 is processed one reference sequence at a time, so whole-genome
DANTE output (gigabytes) runs in bounded memory.

Usage:
    utils/measure_core_seeding.py -g DANTE.gff3
    utils/measure_core_seeding.py -g DANTE.gff3 --sweep_similarity 0.2,0.3,0.4
"""

import argparse
import math
import re
import sys
from collections import Counter, defaultdict

CORE = ('RT', 'RH', 'INT')

# ordered core triplet, in element orientation -> superfamily
SF_BY_ORDER = {
    ('RT', 'RH', 'INT'): 'Ty3_gypsy',
    ('INT', 'RT', 'RH'): 'Ty1_copia',
    }

DB_POS_RE = re.compile(r'^(\d+):(\d+)of(\d+)$')

# how far two annotations may overlap on the genome and still
# be considered fragments of one domain rather than rival hits
MAX_GENOME_OVERLAP = 30


class Domain(object):
    __slots__ = ('start', 'end', 'strand', 'name', 'cls', 'db', 'sim',
                 'ident', 'rel', 'interrupt', 'nfrag')

    def __init__(self, start, end, strand, name, cls, db, sim, ident, rel,
                 interrupt):
        self.start = start
        self.end = end
        self.strand = strand
        self.name = name
        self.cls = cls
        self.db = db
        self.sim = sim
        self.ident = ident
        self.rel = rel
        self.interrupt = interrupt
        self.nfrag = 1


def as_float(value, default=None):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def parse_attributes(field):
    out = {}
    for item in field.rstrip(';').split(';'):
        if not item:
            continue
        key, _, value = item.partition('=')
        out[key.strip()] = value
    return out


def parse_db_pos(value):
    m = DB_POS_RE.match(value or '')
    if not m:
        return None
    return int(m.group(1)), int(m.group(2)), int(m.group(3))


def load_sequences(path):
    """Parse the GFF3 once into per-sequence lists of field tuples.

    Attribute parsing dominates the runtime and the sweeps need many
    passes, so the parse is done once and each pass rebuilds cheap
    Domain objects from the cached tuples.  Repeated strings are
    interned -- domain names and classifications are drawn from a small
    vocabulary, so this is most of the memory back.

    DANTE emits features grouped by reference sequence; a sequence seen
    in more than one block is reported rather than silently split.
    """
    blocks = []
    current = None
    buf = []
    seen = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            col = line.rstrip('\n').split('\t')
            if len(col) < 9 or col[2] != 'protein_domain':
                continue
            seqid = col[0]
            if seqid != current:
                if current is not None:
                    blocks.append((current, buf))
                    if current in seen:
                        sys.stderr.write(
                            'WARNING: {} is not contiguous in the file; '
                            'its domains were processed in separate '
                            'blocks\n'.format(current))
                    seen.add(current)
                current = seqid
                buf = []
            attr = parse_attributes(col[8])
            buf.append((
                int(col[3]), int(col[4]), sys.intern(col[6]),
                sys.intern(attr.get('Name', '')),
                sys.intern(attr.get('Final_Classification', '')),
                attr.get('Best_Hit_DB_Pos', ''),
                as_float(attr.get('Similarity'), 0.0),
                as_float(attr.get('Identity'), 0.0),
                as_float(attr.get('Relat_Length'), 0.0),
                as_float(attr.get('Relat_Interruptions'), 0.0)))
    if current is not None:
        blocks.append((current, buf))
    return blocks


def is_ltr(domain):
    """Domain is an LTR retrotransposon domain by its own classification."""
    return domain.cls.replace('|', '/').startswith('Class_I/LTR')


def passes_filter(domain, args):
    """Primary arm of dante_filtering() (ltr_utils.R:951).

    The secondary arm (`Similarity * Relat_Length > 0.35` for domains
    with a same-cluster neighbour) is not reproduced: it depends on the
    clustering, which is what core mode replaces.  Counting only the
    primary arm makes this a *lower* bound on what core mode would see.
    """
    return (domain.sim >= args.min_similarity and
            domain.ident >= args.min_identity and
            domain.rel >= args.min_relative_length and
            domain.interrupt <= args.max_relat_interruptions)


def merge_split_domains(domains, args, stats):
    """Collapse same-name adjacent annotations that tile one reference
    domain -- the frameshift signature.  `domains` must be sorted by
    coordinate; strand is respected.
    """
    out = []
    for dom in domains:
        prev = out[-1] if out else None
        if (prev is not None and prev.name == dom.name and
                prev.strand == dom.strand):
            gap = dom.start - prev.end - 1
            pa, pb = parse_db_pos(prev.db), parse_db_pos(dom.db)
            # a substantially overlapping pair is two competing hits to
            # the same locus, not one domain broken in two.  The real
            # pipeline drops those in gff_cleanup_overlaps() first.
            if (pa is not None and pb is not None and
                    -MAX_GENOME_OVERLAP <= gap <= args.split_max_gap):
                # on the minus strand the element runs right-to-left, so
                # the reference-order test flips
                if prev.strand == '+':
                    ordered = pb[0] >= pa[0]
                else:
                    ordered = pa[0] >= pb[0]
                overlap = min(pa[1], pb[1]) - max(pa[0], pb[0]) + 1
                shorter = min(pa[1] - pa[0] + 1, pb[1] - pb[0] + 1)
                frac = max(0, overlap) / float(shorter) if shorter else 1.0
                if ordered and frac <= args.split_max_db_overlap:
                    stats['split_pairs'] += 1
                    stats['split_by_name'][dom.name] += 1
                    stats['split_gaps'].append(gap)
                    stats['split_rel_sum'].append(prev.rel + dom.rel)
                    stats['split_min_rel'].append(min(prev.rel, dom.rel))
                    prev.end = max(prev.end, dom.end)
                    prev.rel = min(1.0, prev.rel + dom.rel)
                    prev.sim = max(prev.sim, dom.sim)
                    prev.nfrag += 1
                    continue
                stats['dup_pairs'] += 1
                stats['dup_by_name'][dom.name] += 1
        out.append(dom)
    return out


def find_seeds(domains, args):
    """Enumerate ordered RT/RH/INT seeds (design §6.2).

    `domains` is one reference sequence, sorted by coordinate.  Returns a
    list of (superfamily, indices, span) over the passed-in list.
    """
    seeds = []
    for strand in ('+', '-'):
        idx = [i for i, d in enumerate(domains)
               if d.strand == strand and d.name in CORE and is_ltr(d)]
        if len(idx) < 3:
            continue
        # runs of core candidates separated by at most core_max_gap
        runs = []
        run = [idx[0]]
        for a, b in zip(idx, idx[1:]):
            if domains[b].start - domains[a].end - 1 <= args.core_max_gap:
                run.append(b)
            else:
                runs.append(run)
                run = [b]
        runs.append(run)

        for run in runs:
            if len(run) < 3:
                continue
            order = run if strand == '+' else list(reversed(run))
            cands = []
            n = len(order)
            for i in range(n):
                for j in range(i + 1, n):
                    if _gap(domains, order[i], order[j], strand) > \
                            args.core_max_gap:
                        break
                    for k in range(j + 1, n):
                        if _gap(domains, order[j], order[k], strand) > \
                                args.core_max_gap:
                            break
                        trio = (domains[order[i]].name,
                                domains[order[j]].name,
                                domains[order[k]].name)
                        sf = SF_BY_ORDER.get(trio)
                        if sf is None:
                            continue
                        members = (order[i], order[j], order[k])
                        lo = min(domains[m].start for m in members)
                        hi = max(domains[m].end for m in members)
                        span = hi - lo + 1
                        if span > args.core_max_span:
                            continue
                        # design 6.2 step 3: most compact triple wins,
                        # ties broken on summed Similarity
                        qual = sum(domains[m].sim for m in members)
                        cands.append((span, -qual, members, sf))
            # greedy: most compact first (design 6.2 step 3)
            cands.sort(key=lambda c: (c[0], c[1]))
            used = set()
            for span, _q, members, sf in cands:
                if used.intersection(members):
                    continue
                used.update(members)
                seeds.append((sf, members, span))
    return seeds


def _gap(domains, a, b, strand):
    if strand == '+':
        return domains[b].start - domains[a].end - 1
    return domains[a].start - domains[b].end - 1


def nearest_outside(domains, members, args):
    """Names/distances of the first filtered domain beyond each seed end.

    Answers: would the outward LTR search be blocked immediately?  Under
    the design's blocking-by-default walk (6.3), a leftover core domain
    -- e.g. the second half of a split RT -- stops the window dead.

    `domains` is sorted by start, so the search walks outward from the
    seed's own indices rather than scanning the whole sequence: with
    10^5 domains on a pseudomolecule the full scan is what makes this
    measurement intractable.
    """
    lo = min(domains[m].start for m in members)
    hi = max(domains[m].end for m in members)
    member_set = set(members)

    left = None
    i = min(members) - 1
    while i >= 0:
        if i not in member_set and domains[i].end < lo:
            left = (lo - domains[i].end - 1, domains[i].name)
            break
        i -= 1

    right = None
    j = max(members) + 1
    while j < len(domains):
        if j not in member_set and domains[j].start > hi:
            right = (domains[j].start - hi - 1, domains[j].name)
            break
        j += 1

    return left, right


def quantile(values, q):
    if not values:
        return None
    ordered = sorted(values)
    if len(ordered) == 1:
        return float(ordered[0])
    pos = (len(ordered) - 1) * q
    lo, hi = int(math.floor(pos)), int(math.ceil(pos))
    if lo == hi:
        return float(ordered[lo])
    return ordered[lo] + (ordered[hi] - ordered[lo]) * (pos - lo)


def new_stats():
    return {
        'domains_total': 0,
        'domains_passed': 0,
        'core_passed': Counter(),
        'split_pairs': 0,
        'dup_pairs': 0,
        'split_by_name': Counter(),
        'dup_by_name': Counter(),
        'split_gaps': [],
        'split_rel_sum': [],
        'split_min_rel': [],
        'seeds': Counter(),
        'seed_spans': [],
        'blocked_near': Counter(),
        'blocked_near_core': 0,
        'seeds_total': 0,
        }


def run_pass(blocks, args, merge, dump=None):
    """One evaluation of the cached domains at the current thresholds.

    When `dump` is a list, every seed's location is appended to it --
    used to pick reproducible test fixtures.
    """
    stats = new_stats()
    for seqid, rows in blocks:
        stats['domains_total'] += len(rows)
        kept = [Domain(*row) for row in rows]
        kept.sort(key=lambda d: (d.start, d.end))
        # merge BEFORE filtering: each fragment of a split domain has a
        # reduced Relat_Length by construction, so filtering first
        # destroys exactly the annotations we are trying to rejoin
        if merge:
            kept = merge_split_domains(kept, args, stats)
        kept = [d for d in kept if passes_filter(d, args)]
        stats['domains_passed'] += len(kept)
        for d in kept:
            if d.name in CORE and is_ltr(d):
                stats['core_passed'][d.name] += 1

        seeds = find_seeds(kept, args)
        stats['seeds_total'] += len(seeds)
        if dump is not None:
            for sf, members, span in seeds:
                dump.append((seqid,
                             min(kept[m].start for m in members),
                             max(kept[m].end for m in members),
                             kept[members[0]].strand, sf, span))
        for sf, members, span in seeds:
            stats['seeds'][sf] += 1
            stats['seed_spans'].append(span)
            left, right = nearest_outside(kept, members, args)
            for side in (left, right):
                if side is None:
                    continue
                dist, name = side
                if dist <= args.block_near_bp:
                    stats['blocked_near'][name] += 1
                    if name in CORE:
                        stats['blocked_near_core'] += 1
    return stats


def report(stats, args, label, out=sys.stdout):
    out.write('\n=== {} ===\n'.format(label))
    out.write('domains in file                  : {}\n'
              .format(stats['domains_total']))
    out.write('passing filter                   : {} ({:.1%})\n'
              .format(stats['domains_passed'],
                      stats['domains_passed'] /
                      float(max(stats['domains_total'], 1))))
    out.write('LTR core domains passing         : {}\n'
              .format(', '.join('{} {}'.format(n, stats['core_passed'][n])
                                for n in CORE)))
    if stats['split_pairs'] or stats['dup_pairs']:
        out.write('same-name adjacent pairs         : {} split / '
                  '{} duplicate\n'.format(stats['split_pairs'],
                                          stats['dup_pairs']))
        out.write('  split by domain                : {}\n'.format(
            ', '.join('{} {}'.format(k, v) for k, v in
                      stats['split_by_name'].most_common(8)) or '-'))
        out.write('  duplicate by domain            : {}\n'.format(
            ', '.join('{} {}'.format(k, v) for k, v in
                      stats['dup_by_name'].most_common(8)) or '-'))
        gq = [quantile(stats['split_gaps'], q) for q in (0.5, 0.95, 0.99, 1.0)]
        out.write('  split gap bp (q50 q95 q99 max) : {}\n'.format(
            ' '.join('-' if v is None else '{:.0f}'.format(v) for v in gq)))
        # where the gaps concentrate distinguishes a frameshift (tens of
        # bp) from a domain interrupted by a nested insertion (kb)
        edges = (0, 20, 50, 100, 200, 500, 1000, 2000, 4000)
        hist = [0] * len(edges)
        for g in stats['split_gaps']:
            for bi in range(len(edges) - 1, -1, -1):
                if g >= edges[bi]:
                    hist[bi] += 1
                    break
        out.write('  split gap histogram            : {}\n'.format(
            ', '.join('{}-{}: {}'.format(
                edges[bi],
                edges[bi + 1] if bi + 1 < len(edges) else '+',
                hist[bi]) for bi in range(len(edges)))))
        rq = [quantile(stats['split_rel_sum'], q) for q in (0.05, 0.5, 0.95)]
        mq = [quantile(stats['split_min_rel'], q) for q in (0.05, 0.5, 0.95)]
        out.write('  summed Relat_Length (q05/50/95): {}\n'.format(
            ' '.join('-' if v is None else '{:.2f}'.format(v) for v in rq)))
        out.write('  weaker fragment Relat_Length   : {}\n'.format(
            ' '.join('-' if v is None else '{:.2f}'.format(v) for v in mq)))
    out.write('ordered core seeds               : {} (gypsy {} / copia {})\n'
              .format(stats['seeds_total'], stats['seeds']['Ty3_gypsy'],
                      stats['seeds']['Ty1_copia']))
    sq = [quantile(stats['seed_spans'], q) for q in (0.5, 0.95, 0.99, 1.0)]
    out.write('  seed span (q50 q95 q99 max)    : {}\n'.format(
        ' '.join('-' if v is None else '{:.0f}'.format(v) for v in sq)))
    out.write('seeds whose outward walk meets a domain within {} bp: {}\n'
              .format(args.block_near_bp, sum(stats['blocked_near'].values())))
    out.write('  of which a CORE domain (would block): {}\n'
              .format(stats['blocked_near_core']))
    out.write('  nearest-domain names           : {}\n'.format(
        ', '.join('{} {}'.format(k, v) for k, v in
                  stats['blocked_near'].most_common(8)) or '-'))


def get_arguments():
    parser = argparse.ArgumentParser(
        description=('Measure core-mode seeding on raw DANTE output: '
                     'domain fragmentation, filter thresholds, and how '
                     'many ordered RT/RH/INT seeds result.'),
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-g', '--gff3', required=True,
                        help='raw DANTE GFF3 (not dante_ltr output)')
    parser.add_argument('--min_similarity', type=float, default=0.4)
    parser.add_argument('--min_identity', type=float, default=0.2)
    parser.add_argument('--min_relative_length', type=float, default=0.6)
    parser.add_argument('--max_relat_interruptions', type=float, default=8.0)
    parser.add_argument('--core_max_gap', type=int, default=1500)
    parser.add_argument('--core_max_span', type=int, default=5000)
    parser.add_argument(
        '--split_max_gap', type=int, default=50,
        help=('max bp between two fragments of one domain.  The\n'
              'frameshift excess is concentrated below ~20 bp;\n'
              'complementary tiling stops being specific beyond a\n'
              'few hundred bp (default %(default)s)'))
    parser.add_argument('--split_max_db_overlap', type=float, default=0.3)
    parser.add_argument(
        '--block_near_bp', type=int, default=500,
        help=('a domain this close beyond a seed end truncates the LTR\n'
              'search window almost to nothing (default %(default)s)'))
    parser.add_argument(
        '--dump_seeds', default=None, metavar='TSV',
        help=('write every seed as seqid/start/end/strand/superfamily/span\n'
              'to this file.  Used to choose test fixtures reproducibly.'))
    parser.add_argument(
        '--no_merge', action='store_true',
        help='skip the with/without-merging comparison')
    parser.add_argument(
        '--sweep_similarity', default=None,
        help=('comma-separated Similarity thresholds to sweep, e.g.\n'
              '0.2,0.3,0.4 -- reports seed counts at each'))
    parser.add_argument(
        '--sweep_relative_length', default=None,
        help='comma-separated Relat_Length thresholds to sweep')
    return parser.parse_args()


def main():
    args = get_arguments()

    blocks = load_sequences(args.gff3)
    print('loaded {} domains across {} reference sequences'.format(
        sum(len(rows) for _s, rows in blocks), len(blocks)))

    dump = [] if args.dump_seeds else None
    base = run_pass(blocks, args, merge=False, dump=dump)
    if args.dump_seeds:
        with open(args.dump_seeds, 'w') as fh:
            fh.write('seqid\tstart\tend\tstrand\tsuperfamily\tspan\n')
            for row in dump:
                fh.write('\t'.join(str(x) for x in row) + '\n')
        print('seeds -> {}'.format(args.dump_seeds))
    report(base, args, 'no fragment merging, Similarity >= {} '
                       'Relat_Length >= {}'.format(args.min_similarity,
                                                   args.min_relative_length))

    if not args.no_merge:
        merged = run_pass(blocks, args, merge=True)
        report(merged, args, 'WITH fragment merging, same thresholds')
        delta = merged['seeds_total'] - base['seeds_total']
        print('\nmerging split domains changes seed count by {:+d} '
              '({:+.1%})'.format(
                  delta, delta / float(max(base['seeds_total'], 1))))
        print('seeds blocked at <= {} bp by a core domain: {} -> {}'
              .format(args.block_near_bp, base['blocked_near_core'],
                      merged['blocked_near_core']))

    for flag, attr, tag in (
            (args.sweep_similarity, 'min_similarity', 'Similarity'),
            (args.sweep_relative_length, 'min_relative_length',
             'Relat_Length')):
        if not flag:
            continue
        print('\n{} sweep (with fragment merging)'.format(tag))
        print('{:>8}  {:>10}  {:>8}  {:>8}  {:>8}'.format(
            tag, 'domains', 'seeds', 'gypsy', 'copia'))
        original = getattr(args, attr)
        for raw in flag.split(','):
            value = float(raw)
            setattr(args, attr, value)
            st = run_pass(blocks, args, merge=True)
            print('{:>8.2f}  {:>10}  {:>8}  {:>8}  {:>8}'.format(
                value, st['domains_passed'], st['seeds_total'],
                st['seeds']['Ty3_gypsy'], st['seeds']['Ty1_copia']))
        setattr(args, attr, original)


if __name__ == '__main__':
    main()
