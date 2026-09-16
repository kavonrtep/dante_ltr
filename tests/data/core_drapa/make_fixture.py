#!/usr/bin/env python3
"""Rebuild the core_drapa test fixture from the source Drapa assembly.

The fixture is a 73 kb window of Draparnaldia carrying four ordered
Ty3/gypsy cores, each with a CHD immediately 3' of the core -- the
common layout on this genome (design doc §5.4) and the one that
exercises the traversal rule of §6.3.

Draparnaldia is not covered by REXdb at all: lineage mode finds zero
complete elements on the whole genome, which is what makes this the
fixture that encodes the point of core mode.

The window was chosen by scanning seeds from
`utils/measure_core_seeding.py --dump_seeds` for the most compact region
containing >=2 gypsy seeds that all have a 3' CHD.  Sources live outside
the repository, so this script is for regenerating the fixture, not for
running in CI.

Usage (paths are the defaults recorded in README.md):
    tests/data/core_drapa/make_fixture.py
"""

import argparse
import os

SRC = ('/mnt/ceph/454_data/Drapa/hifiasm/assembly_2025_07_30/'
       'DRA_2025_07_30/output')
DEFAULT_FASTA = SRC + '/hifiasm_assembly.bp.p_ctg.gfa.fasta'
DEFAULT_GFF = SRC + '/analysis/repeat_annotation/output/DANTE/DANTE.gff3'

CONTIG = 'ptg000002l'
START = 19026193          # 1-based, inclusive
END = 19099558
NEW_NAME = 'drapa_ctg1'   # matches the smoke fixture's naming convention


def read_region(fasta, contig, start, end):
    """Pull one 1-based inclusive region using the .fai index."""
    fai = fasta + '.fai'
    offset = line_bases = line_width = length = None
    with open(fai) as fh:
        for line in fh:
            name, ln, off, lb, lw = line.split('\t')[:5]
            if name == contig:
                length, offset = int(ln), int(off)
                line_bases, line_width = int(lb), int(lw)
                break
    if offset is None:
        raise SystemExit('contig {} not found in {}'.format(contig, fai))
    if end > length:
        raise SystemExit('region end {} beyond contig length {}'
                         .format(end, length))

    def seek_of(pos):                      # pos is 0-based
        return offset + pos // line_bases * line_width + pos % line_bases

    with open(fasta) as fh:
        fh.seek(seek_of(start - 1))
        # read generously, then strip newlines and trim to the exact span
        want = end - start + 1
        raw = fh.read(want + want // line_bases + 16)
    return raw.replace('\n', '')[:want]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--fasta', default=DEFAULT_FASTA)
    ap.add_argument('--gff3', default=DEFAULT_GFF)
    ap.add_argument('--outdir', default=os.path.dirname(
        os.path.abspath(__file__)))
    args = ap.parse_args()

    seq = read_region(args.fasta, CONTIG, START, END)
    fasta_out = os.path.join(args.outdir, 'genome.fasta')
    with open(fasta_out, 'w') as fh:
        fh.write('>{}\n'.format(NEW_NAME))
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i + 60] + '\n')

    shift = START - 1
    kept = 0
    gff_out = os.path.join(args.outdir, 'dante.gff3')
    with open(args.gff3) as fin, open(gff_out, 'w') as fout:
        fout.write('##gff-version 3\n')
        for line in fin:
            if line.startswith('#'):
                continue
            col = line.rstrip('\n').split('\t')
            if len(col) < 9 or col[0] != CONTIG:
                continue
            s, e = int(col[3]), int(col[4])
            if s < START or e > END:       # keep only fully contained features
                continue
            col[0] = NEW_NAME
            col[3], col[4] = str(s - shift), str(e - shift)
            fout.write('\t'.join(col) + '\n')
            kept += 1

    print('{}:{}-{} -> {} ({} bp, {} domains)'.format(
        CONTIG, START, END, NEW_NAME, len(seq), kept))


if __name__ == '__main__':
    main()
