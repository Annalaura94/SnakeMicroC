#!/usr/bin/env python3
import argparse
import itertools
import numpy as np
import pandas as pd
from hicrep import hicrepSCC
from hicrep.utils import readMcool


def read_chroms(view_path):
    chroms = []
    with open(view_path) as fh:
        for line in fh:
            if line.strip() and not line.startswith('#'):
                chroms.append(line.split('\t', 1)[0])
    return np.array(chroms, dtype=str)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mcools', nargs='+', required=True)
    ap.add_argument('--samples', nargs='+', required=True)
    ap.add_argument('--view', required=True)
    ap.add_argument('--resolution', type=int, required=True)
    ap.add_argument('--h', type=int, default=1)
    ap.add_argument('--max-distance', type=int, default=5000000)
    ap.add_argument('--out-long', required=True)
    ap.add_argument('--out-matrix', required=True)
    args = ap.parse_args()
    if len(args.mcools) != len(args.samples):
        raise ValueError('mcools and samples must have equal lengths')

    chroms = read_chroms(args.view)
    coolers = {}
    for path, sample in zip(args.mcools, args.samples):
        c, binsize = readMcool(path, args.resolution)
        if int(binsize) != args.resolution:
            raise ValueError(f'{sample}: expected {args.resolution}, got {binsize}')
        coolers[sample] = c

    rows = []
    matrix = pd.DataFrame(np.eye(len(args.samples)), index=args.samples, columns=args.samples)
    for a, b in itertools.combinations(args.samples, 2):
        scc = np.asarray(
            hicrepSCC(
                coolers[a], coolers[b], args.h, args.max_distance,
                False, chroms
            ), dtype=float
        )
        valid = np.isfinite(scc)
        mean_scc = float(np.nanmean(scc)) if valid.any() else np.nan
        matrix.loc[a, b] = mean_scc
        matrix.loc[b, a] = mean_scc
        for chrom, value in zip(chroms, scc):
            rows.append({'sample1': a, 'sample2': b, 'chrom': chrom, 'SCC': value})
        rows.append({'sample1': a, 'sample2': b, 'chrom': 'MEAN', 'SCC': mean_scc})

    pd.DataFrame(rows).to_csv(args.out_long, sep='\t', index=False)
    matrix.to_csv(args.out_matrix, sep='\t')


if __name__ == '__main__':
    main()
