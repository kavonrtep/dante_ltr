#!/usr/bin/env python3
"""Compare lineage-mode and core-mode dante_ltr annotations.

This is the concordance measurement of docs/core_domain_mode_design.md
§10.  Run both modes on the same input and point this at the two GFF3s.

The headline number has to be *conditional on the core being present*.
Core mode requires all three of RT/RH/INT to survive the domain filter,
and on a well-covered genome 2-25% of validated elements do not carry
all three (design §5.2).  Core mode cannot find those, by construction,
so an unconditional recovery rate measures the filter, not the
algorithm.  Both figures are reported.

Usage:
    utils/compare_detection_modes.py -l lineage.gff3 -c core.gff3
    utils/compare_detection_modes.py -l lineage.gff3 -c core.gff3 \\
        --tolerance 20 --dante DANTE.gff3
"""

import argparse
import sys
from collections import Counter, defaultdict

CORE_DOMAINS = ('RT', 'RH', 'INT')
RANKS = ('D', 'DL', 'DLT', 'DLP', 'DLTP')
COMPLETE_RANKS = ('DL', 'DLT', 'DLP', 'DLTP')


def parse_attributes(field):
    out = {}
    for item in field.rstrip(';').split(';'):
        if not item:
            continue
        key, _, value = item.partition('=')
        out[key.strip()] = value
    return out


def superfamily_of(cls):
    if not cls:
        return None
    parts = cls.replace('|', '/').split('/')
    if len(parts) < 3 or parts[0] != 'Class_I' or parts[1] != 'LTR':
        return None
    if parts[2].startswith('Ty1'):
        return 'Ty1/copia'
    if parts[2].startswith('Ty3'):
        return 'Ty3/gypsy'
    return None


def depth(cls):
    return len(cls.split('|')) if cls else 0


class Element(object):
    __slots__ = ('seqid', 'start', 'end', 'strand', 'rank', 'cls', 'eid',
                 'domains')

    def __init__(self, seqid, start, end, strand, attr):
        self.seqid = seqid
        self.start = start
        self.end = end
        self.strand = strand
        self.rank = attr.get('Rank', '')
        self.cls = attr.get('Final_Classification', '')
        self.eid = attr.get('ID', '')
        self.domains = []

    def key(self):
        return (self.seqid, self.strand)


def read_elements(path):
    """Elements plus the domain names each one contains."""
    elements = {}
    children = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            col = line.rstrip('\n').split('\t')
            if len(col) < 9:
                continue
            attr = parse_attributes(col[8])
            if col[2] == 'transposable_element':
                eid = attr.get('ID')
                if eid:
                    elements[eid] = Element(col[0], int(col[3]), int(col[4]),
                                            col[6], attr)
            elif col[2] == 'protein_domain':
                parent = attr.get('Parent')
                if parent:
                    children[parent].append(attr.get('Name', ''))
    for eid, names in children.items():
        if eid in elements:
            elements[eid].domains = names
    return list(elements.values())


def has_complete_core(el):
    return set(CORE_DOMAINS).issubset(set(el.domains))


def match_elements(a_list, b_list, tolerance):
    """Greedy nearest-boundary matching within `tolerance` bp.

    Returns (pairs, a_only, b_only).
    """
    by_key = defaultdict(list)
    for j, b in enumerate(b_list):
        by_key[b.key()].append(j)
    used = set()
    pairs = []
    a_only = []
    for a in a_list:
        best, best_d = None, None
        for j in by_key.get(a.key(), ()):
            if j in used:
                continue
            b = b_list[j]
            d = abs(a.start - b.start) + abs(a.end - b.end)
            if d <= 2 * tolerance and (best_d is None or d < best_d):
                best, best_d = j, d
        if best is None:
            a_only.append(a)
        else:
            used.add(best)
            pairs.append((a, b_list[best], best_d))
    b_only = [b for j, b in enumerate(b_list) if j not in used]
    return pairs, a_only, b_only


def pct(n, d):
    return '{:.1f}%'.format(100.0 * n / d) if d else '-'


def rank_table(name, els, out):
    c = Counter(e.rank for e in els)
    out.write('  {:<10}'.format(name))
    for r in RANKS:
        out.write('{:>8}'.format(c.get(r, 0)))
    out.write('{:>10}\n'.format(len(els)))


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-l', '--lineage', required=True,
                    help='dante_ltr output, --mode lineage')
    ap.add_argument('-c', '--core', required=True,
                    help='dante_ltr output, --mode core')
    ap.add_argument('--tolerance', type=int, default=20,
                    help='per-boundary tolerance in bp (default %(default)s)')
    args = ap.parse_args()

    lin = read_elements(args.lineage)
    cor = read_elements(args.core)
    lin_c = [e for e in lin if e.rank in COMPLETE_RANKS]
    cor_c = [e for e in cor if e.rank in COMPLETE_RANKS]

    out = sys.stdout
    out.write('\n=== elements by rank ===\n')
    out.write('  {:<10}'.format(''))
    for r in RANKS:
        out.write('{:>8}'.format(r))
    out.write('{:>10}\n'.format('total'))
    rank_table('lineage', lin, out)
    rank_table('core', cor, out)

    out.write('\n=== concordance (complete elements, +/-{} bp) ===\n'
              .format(args.tolerance))
    pairs, lin_only, cor_only = match_elements(lin_c, cor_c, args.tolerance)
    out.write('  lineage complete elements      : {}\n'.format(len(lin_c)))
    out.write('  core    complete elements      : {}\n'.format(len(cor_c)))
    out.write('  matched                        : {}  ({} of lineage)\n'
              .format(len(pairs), pct(len(pairs), len(lin_c))))
    out.write('  lineage-only                   : {}\n'.format(len(lin_only)))
    out.write('  core-only                      : {}\n'.format(len(cor_only)))

    # The conditional figure: core mode can only find elements whose
    # three core domains survived the filter.
    eligible = [e for e in lin_c if has_complete_core(e)]
    matched_eligible = sum(1 for a, _b, _d in pairs if has_complete_core(a))
    out.write('\n  lineage elements WITH a complete filtered core: {} ({} of'
              ' all complete)\n'.format(len(eligible),
                                        pct(len(eligible), len(lin_c))))
    out.write('  of those, recovered by core mode             : {}  ({})'
              '   <-- the design §10 gate\n'
              .format(matched_eligible, pct(matched_eligible, len(eligible))))

    missed_no_core = sum(1 for e in lin_only if not has_complete_core(e))
    out.write('  lineage-only explained by an incomplete core : {} of {}\n'
              .format(missed_no_core, len(lin_only)))

    if pairs:
        exact = sum(1 for _a, _b, d in pairs if d == 0)
        out.write('\n  boundary agreement on matched pairs:\n')
        out.write('    exact (0 bp)                 : {}  ({})\n'
                  .format(exact, pct(exact, len(pairs))))
        for lim in (2, 10, 20):
            n = sum(1 for _a, _b, d in pairs if d <= 2 * lim)
            out.write('    within +/-{:<3} bp per end      : {}  ({})\n'
                      .format(lim, n, pct(n, len(pairs))))

    out.write('\n=== superfamily agreement on matched pairs ===\n')
    agree = disagree = unknown = 0
    for a, b, _d in pairs:
        sa, sb = superfamily_of(a.cls), superfamily_of(b.cls)
        if sa is None or sb is None:
            unknown += 1
        elif sa == sb:
            agree += 1
        else:
            disagree += 1
            if disagree <= 5:
                out.write('    MISMATCH {}:{}-{}  lineage={}  core={}\n'
                          .format(a.seqid, a.start, a.end, a.cls, b.cls))
    out.write('  agree {}   disagree {}   undetermined {}\n'
              .format(agree, disagree, unknown))
    if disagree == 0:
        out.write('  the order-derived superfamily never contradicts '
                  'lineage mode\n')

    out.write('\n=== classification depth on matched pairs ===\n')
    same = deeper = shallower = 0
    for a, b, _d in pairs:
        da, db = depth(a.cls), depth(b.cls)
        if db == da:
            same += 1
        elif db > da:
            deeper += 1
        else:
            shallower += 1
    out.write('  core same depth as lineage     : {}  ({})\n'
              .format(same, pct(same, len(pairs))))
    out.write('  core shallower (demoted)       : {}  ({})\n'
              .format(shallower, pct(shallower, len(pairs))))
    out.write('  core deeper                    : {}\n'.format(deeper))

    out.write('\n=== core-only elements ===\n')
    out.write('  count                          : {}\n'.format(len(cor_only)))
    if cor_only:
        c = Counter(superfamily_of(e.cls) or 'unclassified' for e in cor_only)
        out.write('  by superfamily                 : {}\n'.format(
            ', '.join('{} {}'.format(k, v) for k, v in c.most_common())))
        c = Counter(e.rank for e in cor_only)
        out.write('  by rank                        : {}\n'.format(
            ', '.join('{} {}'.format(k, c[k]) for k in RANKS if k in c)))
    out.write('\n')


if __name__ == '__main__':
    main()
