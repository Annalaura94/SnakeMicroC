#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def parse_stats(path):
    stats = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            key = parts[0]
            try:
                value = float(parts[1])
            except ValueError:
                continue
            stats[key] = value
    return stats


def get(stats, key):
    if key not in stats:
        raise KeyError(f'Missing expected pairtools statistic: {key}')
    return stats[key]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--shallow-qc', action='store_true',
                    help='Pilot screening only: do not apply deep-sequencing yield PASS/FAIL')
    ap.add_argument('--stats', nargs='+', required=True)
    ap.add_argument('--samples', nargs='+', required=True)
    ap.add_argument('--conditions', nargs='+', required=True)
    ap.add_argument('--min-nodup', type=float, required=True)
    ap.add_argument('--min-cis1kb-frac', type=float, required=True)
    ap.add_argument('--recommended-total', type=float, required=True)
    ap.add_argument('--out-tsv', required=True)
    ap.add_argument('--out-pdf', required=True)
    ap.add_argument('--out-status', required=True)
    args = ap.parse_args()

    if not (len(args.stats) == len(args.samples) == len(args.conditions)):
        raise ValueError('stats, samples and conditions must have equal lengths')

    rows = []
    for path, sample, condition in zip(args.stats, args.samples, args.conditions):
        s = parse_stats(path)
        total = get(s, 'total')
        unmapped = get(s, 'total_unmapped')
        mapped = get(s, 'total_mapped')
        dups = get(s, 'total_dups')
        nodup = get(s, 'total_nodups')
        cis = get(s, 'cis')
        trans = get(s, 'trans')
        cis1 = get(s, 'cis_1kb+')
        cis10 = get(s, 'cis_10kb+')
        low_mapq = max(total - unmapped - mapped, 0)
        valid = cis1 + trans
        cis_lt1 = max(cis - cis1, 0)

        cis1_frac = cis1 / nodup if nodup else np.nan
        pass_nodup = nodup > args.min_nodup
        pass_cis1 = cis1_frac > args.min_cis1kb_frac
        qc_pass = bool(pass_nodup and pass_cis1)
        shallow_screen = "SCREEN_PASS" if pass_cis1 else "REVIEW_CIS_FRACTION"

        rows.append({
            'sample': sample,
            'condition': condition,
            'total_read_pairs': total,
            'unmapped_read_pairs': unmapped,
            'low_mapq_read_pairs': low_mapq,
            'mapped_read_pairs': mapped,
            'pcr_dup_read_pairs': dups,
            'no_dup_read_pairs': nodup,
            'no_dup_cis_read_pairs': cis,
            'no_dup_trans_read_pairs': trans,
            'no_dup_valid_read_pairs': valid,
            'no_dup_cis_lt_1kb': cis_lt1,
            'no_dup_cis_ge_1kb': cis1,
            'no_dup_cis_ge_10kb': cis10,
            'mapped_pct_total': 100 * mapped / total if total else np.nan,
            'duplicates_pct_total': 100 * dups / total if total else np.nan,
            'duplicates_pct_mapped': 100 * dups / mapped if mapped else np.nan,
            'nodup_pct_total': 100 * nodup / total if total else np.nan,
            'trans_pct_nodup': 100 * trans / nodup if nodup else np.nan,
            'cis_ge_1kb_pct_nodup': 100 * cis1 / nodup if nodup else np.nan,
            'cis_ge_10kb_pct_nodup': 100 * cis10 / nodup if nodup else np.nan,
            'valid_pct_nodup': 100 * valid / nodup if nodup else np.nan,
            'recommended_total_depth': 'NOT_EVALUATED' if args.shallow_qc else total >= args.recommended_total,
            'pass_nodup_pairs': 'NOT_EVALUATED' if args.shallow_qc else pass_nodup,
            'pass_cis_ge_1kb_fraction': pass_cis1,
            'dovetail_deep_qc_pass': 'NOT_EVALUATED' if args.shallow_qc else qc_pass,
            'qc_mode': 'shallow_qc' if args.shallow_qc else 'full',
            'shallow_screening_status': shallow_screen if args.shallow_qc else 'NOT_APPLICABLE',
        })

    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(args.out_tsv), exist_ok=True)
    df.to_csv(args.out_tsv, sep='\t', index=False)

    with open(args.out_status, 'w') as out:
        for _, r in df.iterrows():
            out.write(
                f"{r['sample']}\t{r['shallow_screening_status'] if args.shallow_qc else ('PASS' if r['dovetail_deep_qc_pass'] else 'FAIL')}"
                f"\tNoDup={int(r['no_dup_read_pairs'])}"
                f"\tCis>=1kb={r['cis_ge_1kb_pct_nodup']:.1f}%\n"
            )

    order = list(df['sample'])
    x = np.arange(len(order))

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))

    axes[0].bar(x, df['mapped_pct_total'], label='Mapped')
    axes[0].bar(x, df['low_mapq_read_pairs'] / df['total_read_pairs'] * 100,
                bottom=df['mapped_pct_total'], label='Low MAPQ')
    axes[0].bar(x, df['unmapped_read_pairs'] / df['total_read_pairs'] * 100,
                bottom=df['mapped_pct_total'] + df['low_mapq_read_pairs'] / df['total_read_pairs'] * 100,
                label='Unmapped')
    axes[0].set_title('Alignment structure')
    axes[0].set_ylabel('% total read pairs')
    axes[0].set_ylim(0, 100)

    axes[1].bar(x, df['cis_ge_1kb_pct_nodup'], label='cis >=1 kb')
    axes[1].bar(x, df['trans_pct_nodup'], bottom=df['cis_ge_1kb_pct_nodup'], label='trans')
    short = 100 - df['cis_ge_1kb_pct_nodup'] - df['trans_pct_nodup']
    axes[1].bar(x, short, bottom=df['cis_ge_1kb_pct_nodup'] + df['trans_pct_nodup'], label='cis <1 kb')
    axes[1].axhline(args.min_cis1kb_frac * 100, linestyle='--', linewidth=1)
    axes[1].set_title('No-dup contact composition')
    axes[1].set_ylabel('% no-dup read pairs')
    axes[1].set_ylim(0, 100)

    nodup_m = df['no_dup_read_pairs'] / 1e6
    axes[2].bar(x, nodup_m)
    if not args.shallow_qc:
        axes[2].axhline(args.min_nodup / 1e6, linestyle='--', linewidth=1)
    axes[2].set_title('Observed shallow yield' if args.shallow_qc else 'Deep-sequencing yield')
    axes[2].set_ylabel('No-dup read pairs (millions)')

    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(order, rotation=35, ha='right')
        ax.spines[['top', 'right']].set_visible(False)
    axes[0].legend(fontsize=8)
    axes[1].legend(fontsize=8)
    fig.suptitle('Micro-C - shallow pilot QC (depth criterion not evaluated)' if args.shallow_qc else 'Micro-C - QC summary')
    fig.tight_layout()
    fig.savefig(args.out_pdf)
    plt.close(fig)


if __name__ == '__main__':
    main()
