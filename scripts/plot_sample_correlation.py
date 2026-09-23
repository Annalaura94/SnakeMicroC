#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cooler


def vectorize(uri):
    clr = cooler.Cooler(uri)
    mat = np.asarray(clr.matrix(balance=True)[:], dtype=float)
    mat = np.nan_to_num(mat, nan=0.0, posinf=0.0, neginf=0.0)
    iu = np.triu_indices_from(mat, k=1)
    return np.log10(mat[iu] + 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cools', nargs='+', required=True)
    ap.add_argument('--samples', nargs='+', required=True)
    ap.add_argument('--out-tsv', required=True)
    ap.add_argument('--out-pdf', required=True)
    args = ap.parse_args()
    if len(args.cools) != len(args.samples):
        raise ValueError('cools and samples must have equal lengths')

    vectors = [vectorize(x) for x in args.cools]
    lengths = {len(x) for x in vectors}
    if len(lengths) != 1:
        raise ValueError('Contact maps have different bin dimensions')
    df = pd.DataFrame(np.vstack(vectors).T, columns=args.samples)
    corr = df.corr(method='pearson')
    corr.to_csv(args.out_tsv, sep='\t')

    fig, ax = plt.subplots(figsize=(5.5, 5))
    im = ax.imshow(corr.values, vmin=0, vmax=1, origin='lower', aspect='equal')
    ax.set_xticks(range(len(args.samples)))
    ax.set_yticks(range(len(args.samples)))
    ax.set_xticklabels(args.samples, rotation=40, ha='right')
    ax.set_yticklabels(args.samples)
    for i in range(len(args.samples)):
        for j in range(len(args.samples)):
            ax.text(j, i, f'{corr.iloc[i, j]:.3f}', ha='center', va='center', fontsize=8)
    ax.set_title('1 Mb contact-map correlation')
    fig.colorbar(im, ax=ax, label='Pearson r')
    fig.tight_layout()
    fig.savefig(args.out_pdf)
    plt.close(fig)


if __name__ == '__main__':
    main()
