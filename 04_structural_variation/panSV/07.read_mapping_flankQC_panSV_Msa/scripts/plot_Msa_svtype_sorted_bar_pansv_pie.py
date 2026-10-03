#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

SVTYPES = ['DEL', 'INS', 'TRA', 'DUP', 'INV']
LABELS = {
    'DEL': 'Deletion',
    'INS': 'Insertion',
    'TRA': 'Translocation',
    'DUP': 'Duplication',
    'INV': 'Inversion',
}
COLORS = {
    'DEL': '#56B4C3',
    'INS': '#45A68C',
    'TRA': '#5A6682',
    'DUP': '#D88A56',
    'INV': '#D9635F',
}

def read_pansv_type_counts(stats_path: Path):
    stats = pd.read_csv(stats_path, sep='\t', dtype={'section': str, 'key': str, 'count': str})
    rows = stats[stats['section'] == 'SVTYPE'].copy()
    rows['count'] = rows['count'].astype(int)
    counts = {row['key']: int(row['count']) for _, row in rows.iterrows()}
    return {sv: counts.get(sv, 0) for sv in SVTYPES}

def format_k(x, _pos):
    if x >= 1000:
        return f'{int(x/1000)}k'
    return str(int(x))

def autopct_counts(values):
    total = sum(values)
    def _fmt(pct):
        value = int(round(pct * total / 100.0))
        if value == 0:
            return ''
        return f'{value:,}'
    return _fmt

def main():
    ap = argparse.ArgumentParser(description='Draw sorted species-level SV stacked bar and non-redundant panSV pie chart.')
    ap.add_argument('--overview', required=True, type=Path)
    ap.add_argument('--pansv-stats', required=True, type=Path)
    ap.add_argument('--out-prefix', required=True, type=Path)
    args = ap.parse_args()

    df = pd.read_csv(args.overview, sep='\t')
    for sv in SVTYPES:
        if sv not in df.columns:
            df[sv] = 0
    df['total_for_plot'] = df[SVTYPES].sum(axis=1)
    df = df.sort_values(['total_for_plot', 'species'], ascending=[True, True]).reset_index(drop=True)

    pie_counts = read_pansv_type_counts(args.pansv_stats)
    pie_values = [pie_counts[sv] for sv in SVTYPES]

    plt.rcParams.update({
        'font.family': 'DejaVu Sans',
        'font.size': 8.5,
        'axes.linewidth': 1.0,
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
        'svg.fonttype': 'none',
    })

    fig = plt.figure(figsize=(7.2, 4.6))
    ax = fig.add_axes([0.09, 0.16, 0.84, 0.74])

    x = list(range(len(df)))
    bottom = [0] * len(df)
    for sv in SVTYPES:
        vals = df[sv].astype(int).tolist()
        ax.bar(x, vals, bottom=bottom, width=0.62, color=COLORS[sv], edgecolor='none', label=LABELS[sv])
        bottom = [b + v for b, v in zip(bottom, vals)]

    ax.set_xlim(-0.65, len(df) - 0.35)
    ymax = max(bottom) * 1.18
    ax.set_ylim(0, ymax)
    ax.set_ylabel('SV number', fontsize=10)
    ax.set_xlabel('Species', fontsize=9)
    ax.yaxis.set_major_formatter(FuncFormatter(format_k))
    ax.set_xticks(x)
    ax.set_xticklabels(df['species'].tolist(), rotation=45, ha='right', fontsize=7.5)
    ax.tick_params(axis='both', width=1.0, length=3, pad=2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', linestyle='-', alpha=0.12, linewidth=0.6)
    ax.legend(frameon=False, loc='upper left', bbox_to_anchor=(0.01, 0.98), fontsize=8, handlelength=1.1, handletextpad=0.4, labelspacing=0.35)

    ax.text(-0.12, 1.03, 'A', transform=ax.transAxes, fontsize=13, fontweight='bold', va='top', ha='left')

    pie_ax = inset_axes(ax, width='35%', height='47%', loc='upper center', bbox_to_anchor=(0.06, 0.03, 0.88, 0.94), bbox_transform=ax.transAxes, borderpad=0)
    wedges, texts, autotexts = pie_ax.pie(
        pie_values,
        colors=[COLORS[sv] for sv in SVTYPES],
        startangle=90,
        counterclock=False,
        autopct=autopct_counts(pie_values),
        pctdistance=1.20,
        labeldistance=1.05,
        wedgeprops={'linewidth': 0.6, 'edgecolor': 'white'},
        textprops={'fontsize': 7.5, 'color': '#333333'},
    )
    pie_ax.set_aspect('equal')
    pie_ax.set_title('Non-redundant panSV', fontsize=8.5, pad=1.5)

    # Add a small source note, unobtrusive but keeps the figure auditable.
    note = 'Bars: per-species merged SVs after read-mapping flank QC; Pie: non-redundant panSV events.'
    fig.text(0.09, 0.035, note, fontsize=6.8, color='#555555')

    out = args.out_prefix
    out.parent.mkdir(parents=True, exist_ok=True)
    for ext, dpi in [('pdf', None), ('svg', None), ('png', 450)]:
        path = out.with_suffix('.' + ext)
        fig.savefig(path, dpi=dpi, bbox_inches='tight')
    sorted_tsv = out.with_suffix('.sorted_species_counts.tsv')
    df[['species', 'total_for_plot'] + SVTYPES].rename(columns={'total_for_plot': 'total'}).to_csv(sorted_tsv, sep='\t', index=False)
    print(f'WROTE\t{out.with_suffix(".pdf")}')
    print(f'WROTE\t{out.with_suffix(".svg")}')
    print(f'WROTE\t{out.with_suffix(".png")}')
    print(f'WROTE\t{sorted_tsv}')
    print('PIE_COUNTS\t' + '\t'.join(f'{sv}={pie_counts[sv]}' for sv in SVTYPES))

if __name__ == '__main__':
    main()