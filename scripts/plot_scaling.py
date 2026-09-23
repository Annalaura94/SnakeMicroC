#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def load_scaling(path):
    df = pd.read_csv(path, sep='\t', comment='#')
    required = {'min_dist', 'max_dist', 'n_pairs', 'n_bp2'}
    if not required.issubset(df.columns):
        raise ValueError(
            f'{path}: expected pairtools scaling columns {sorted(required)}, '
            f'found {list(df.columns)}'
        )
    for c in required:
        df[c] = pd.to_numeric(df[c], errors='coerce')
    df = df.dropna(subset=list(required)).copy()
    df = df[(df['min_dist'] > 0) & (df['max_dist'] > 0) & (df['n_bp2'] > 0)]
    agg = (
        df.groupby(['min_dist', 'max_dist'], as_index=False)
          .agg(n_pairs=('n_pairs', 'sum'), n_bp2=('n_bp2', 'sum'))
    )
    total = agg['n_pairs'].sum()
    agg['distance_bp'] = np.sqrt(agg['min_dist'] * agg['max_dist'])
    agg['frequency'] = agg['n_pairs'] / agg['n_bp2']
    if total > 0:
        agg['frequency'] = agg['frequency'] / total
    return agg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--inputs', nargs='+', required=True)
    ap.add_argument('--samples', nargs='+', required=True)
    ap.add_argument('--out-pdf', required=True)
    ap.add_argument('--out-tsv', required=True)
    args = ap.parse_args()
    if len(args.inputs) != len(args.samples):
        raise ValueError('inputs and samples must have equal lengths')

    all_df = []
    fig, ax = plt.subplots(figsize=(7, 5))
    for path, sample in zip(args.inputs, args.samples):
        df = load_scaling(path)
        df['sample'] = sample
        all_df.append(df)
        mask = df['frequency'] > 0
        ax.plot(df.loc[mask, 'distance_bp'], df.loc[mask, 'frequency'], label=sample)

    out = pd.concat(all_df, ignore_index=True)
    out.to_csv(args.out_tsv, sep='\t', index=False)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Genomic distance (bp)')
    ax.set_ylabel('Relative cis contact frequency')
    ax.set_title('Micro-C contact decay P(s)')
    ax.legend(fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    fig.savefig(args.out_pdf)
    plt.close(fig)


if __name__ == '__main__':
    main()
