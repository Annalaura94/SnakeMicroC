#!/usr/bin/env python3
import argparse
import re


def canonical_rank(chrom):
    c = chrom[3:] if chrom.startswith('chr') else chrom
    if c.isdigit() and 1 <= int(c) <= 22:
        return (0, int(c))
    if c == 'X':
        return (1, 23)
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fai', required=True)
    ap.add_argument('--chrom-sizes', required=True)
    ap.add_argument('--main-view', required=True)
    ap.add_argument('--autosomes-view', required=True)
    args = ap.parse_args()

    records = []
    with open(args.fai) as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            records.append((f[0], int(f[1])))

    with open(args.chrom_sizes, 'w') as out:
        for chrom, size in records:
            out.write(f'{chrom}\t{size}\n')

    main = [(c, s) for c, s in records if canonical_rank(c) is not None]
    main.sort(key=lambda x: canonical_rank(x[0]))
    autosomes = []
    for c, s in main:
        label = c[3:] if c.startswith('chr') else c
        if label.isdigit() and 1 <= int(label) <= 22:
            autosomes.append((c, s))

    if len(autosomes) != 22:
        raise RuntimeError(
            f'Expected 22 autosomes in FAI, found {len(autosomes)}. '
            'Check chromosome naming in the configured reference.'
        )

    with open(args.main_view, 'w') as out:
        for chrom, size in main:
            out.write(f'{chrom}\t0\t{size}\t{chrom}\n')

    with open(args.autosomes_view, 'w') as out:
        for chrom, size in autosomes:
            out.write(f'{chrom}\t0\t{size}\t{chrom}\n')


if __name__ == '__main__':
    main()
