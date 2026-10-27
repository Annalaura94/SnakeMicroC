#!/usr/bin/env python3
import argparse
import math
import numpy as np
import matplotlib.pyplot as plt
import cooler


def load_matrix(uri, balanced=True):
    clr = cooler.Cooler(uri)
    m = clr.matrix(balance=balanced)[:]
    m = np.asarray(m, dtype=float)
    m = np.nan_to_num(m, nan=0.0, posinf=0.0, neginf=0.0)
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cools', nargs='+', required=True)
    ap.add_argument('--samples', nargs='+', required=True)
    ap.add_argument('--out-pdf', required=True)
    ap.add_argument('--unbalanced', action='store_true')
    ap.add_argument('--clip-percentile', type=float, default=99.7)
    args = ap.parse_args()
    if len(args.cools) != len(args.samples):
        raise ValueError('cools and samples must have equal lengths')

    mats = [load_matrix(x, balanced=not args.unbalanced) for x in args.cools]
    vals = np.concatenate([np.log10(m[m > 0] + 1) for m in mats if np.any(m > 0)])
    vmax = np.percentile(vals, args.clip_percentile) if vals.size else 1.0

    n = len(mats)
    ncols = 2
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(11.5, 5.2 * nrows), constrained_layout=True)
    axes = np.atleast_1d(axes).ravel()
    for ax, mat, sample in zip(axes, mats, args.samples):
        vis = np.log10(mat + 1)
        masked = np.ma.masked_where(mat <= 0, vis)
        im = ax.imshow(masked, origin='lower', interpolation='nearest', aspect='auto', vmin=0, vmax=vmax)
        ax.set_title(sample)
        ax.set_xlabel('1 Mb bins')
        ax.set_ylabel('1 Mb bins')
    for ax in axes[n:]:
        ax.axis('off')
    fig.colorbar(im, ax=axes[:n].tolist(), shrink=0.7, label='log10(contact + 1)')
    fig.suptitle('Micro-C QC maps (raw counts, shallow pilot)' if args.unbalanced else 'Genome-wide Micro-C contact maps - QC view')
    fig.savefig(args.out_pdf, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    main()
