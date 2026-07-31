"""
NIHSS vs. Network FC — within-cluster Spearman across 3 timepoints
-------------------------------------------------------------------
Replicates the acute-session within-cluster analysis at 3 months and 1 year.
Produces:
  - Combined figure: 3 rows (timepoints) × 7 cols (networks)
  - Summary CSV with all rho/p values
"""

import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
import os, warnings
warnings.filterwarnings('ignore')

DATA_PATH = '/home/cemal/Desktop/Opus/stroke_clinical_and_neural_data_with_lesion_disconnectome.pkl'
MAP_PATH  = '/home/cemal/Desktop/Opus/R4_analyses/schaefer400_network_mapping.csv'
OUT_DIR   = '/home/cemal/Desktop/Opus/NIHSS_FC_correlation'

print('Loading data...')
with open(DATA_PATH, 'rb') as f:
    data = pickle.load(f)

net_map = pd.read_csv(MAP_PATH)
networks = ['Vis', 'SomMot', 'DorsAttn', 'SalVentAttn', 'Limbic', 'Cont', 'Default']
net_labels = {
    'Vis':         'Visual',
    'SomMot':      'Somatomotor',
    'DorsAttn':    'Dorsal Attention',
    'SalVentAttn': 'Salience / Ventral Attention',
    'Limbic':      'Limbic',
    'Cont':        'Frontoparietal Control',
    'Default':     'Default Mode',
}
net_to_idx = {net: net_map[net_map['network'] == net]['roi_idx'].values - 1
              for net in networks}

stroke = data[data['subj_type'] == 0].copy()
excluded = stroke.loc[stroke['inclusion_notes'].notna(), 'ID'].unique()
stroke = stroke[~stroke['ID'].isin(excluded)].copy()

sessions = [
    ('acute',     '2 weeks'),
    ('followup',  '3 months'),
    ('followup2', '1 year'),
]

cat_bins   = [-0.1, 1, 5, 10, 100]
cat_labels = ['0–1', '2–5', '5–10', '10+']
cat_colors = {'0–1': '#2ecc71', '2–5': '#f39c12', '5–10': '#e74c3c', '10+': '#8e44ad'}

def within_network_mean(mat, idx):
    sub = mat[np.ix_(idx, idx)]
    return float(np.nanmean(sub[np.triu_indices(len(idx), k=1)]))

# Build per-timepoint DataFrames
session_dfs = {}
for sess_key, sess_label in sessions:
    sess = stroke[stroke['Session'] == sess_key].copy()
    rows = []
    for _, row in sess.iterrows():
        mat   = row['mean_corrmats']
        nihss = row['nih_total']
        if not isinstance(mat, np.ndarray) or mat.shape != (400, 400):
            continue
        if pd.isna(nihss):
            continue
        rec = {'ID': row['ID'], 'nih_total': nihss}
        for net in networks:
            rec[f'fc_{net}'] = within_network_mean(mat, net_to_idx[net])
        rows.append(rec)
    df = pd.DataFrame(rows)
    df['nihss_cat'] = pd.cut(df['nih_total'], bins=cat_bins, labels=cat_labels)
    session_dfs[sess_label] = df
    print(f'\n{sess_label} — n={len(df)}, NIHSS range: {df["nih_total"].min():.0f}–{df["nih_total"].max():.0f}')
    for cat in cat_labels:
        print(f'  {cat}: n={(df["nihss_cat"]==cat).sum()}')

# Combined figure: rows=timepoints, cols=networks
n_tp  = len(sessions)
n_net = len(networks)
fig, axes = plt.subplots(n_tp, n_net, figsize=(32, 14), sharey=False)

all_results = []

for ri, (sess_key, sess_label) in enumerate(sessions):
    df = session_dfs[sess_label]

    for ci, net in enumerate(networks):
        ax     = axes[ri, ci]
        fc_col = f'fc_{net}'
        plot_df = df[['nih_total', fc_col, 'nihss_cat']].dropna()

        # Scatter — colored by cluster
        for cat in cat_labels:
            mask = plot_df['nihss_cat'] == cat
            ax.scatter(
                plot_df.loc[mask, 'nih_total'],
                plot_df.loc[mask, fc_col],
                color=cat_colors[cat], s=40, alpha=0.80,
                edgecolors='white', linewidths=0.3, zorder=3
            )

        # Within-cluster regression + Spearman
        stat_lines = []
        for cat in cat_labels:
            sub = plot_df[plot_df['nihss_cat'] == cat]
            x_c = sub['nih_total'].values
            y_c = sub[fc_col].values
            n_c = len(x_c)

            if n_c >= 3 and x_c.std() > 0:
                rho, pval = stats.spearmanr(x_c, y_c)
                m, b = np.polyfit(x_c, y_c, 1)
                xline = np.linspace(x_c.min(), x_c.max(), 100)
                ax.plot(xline, m * xline + b,
                        color=cat_colors[cat], lw=1.8, ls='--', zorder=4, alpha=0.9)
                stars = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else 'n.s.'
                p_str = 'p<.001' if pval < 0.001 else f'p={pval:.3f}'
                stat_lines.append(f'{cat}: ρ={rho:.2f} {stars} (n={n_c})')
                all_results.append({
                    'timepoint': sess_label, 'network': net_labels[net],
                    'cluster': cat, 'n': n_c,
                    'rho': round(rho, 4), 'p': round(pval, 4), 'sig': stars,
                })
            else:
                stat_lines.append(f'{cat}: n={n_c} (skip)')
                all_results.append({
                    'timepoint': sess_label, 'network': net_labels[net],
                    'cluster': cat, 'n': n_c,
                    'rho': np.nan, 'p': np.nan, 'sig': '-',
                })

        ax.text(0.98, 0.98, '\n'.join(stat_lines),
                transform=ax.transAxes, ha='right', va='top', fontsize=6.5,
                fontfamily='monospace',
                bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='grey', alpha=0.85))

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(labelsize=7)

        if ri == 0:
            ax.set_title(net_labels[net], fontsize=9, fontweight='bold', pad=4)
        if ci == 0:
            ax.set_ylabel(f'{sess_label}\nMean within-network FC (r)', fontsize=8)
        if ri == n_tp - 1:
            ax.set_xlabel('NIHSS score', fontsize=8)

# Shared legend
patches = [mpatches.Patch(color=cat_colors[c], label=c) for c in cat_labels]
fig.legend(handles=patches, title='NIHSS cluster', loc='lower right',
           bbox_to_anchor=(0.99, 0.01), fontsize=9, title_fontsize=10,
           frameon=True, edgecolor='grey')

fig.suptitle(
    'Mean within-network FC vs. NIHSS score — within-cluster Spearman ρ\n'
    '2 weeks  |  3 months  |  1 year post-stroke',
    fontsize=13, fontweight='bold', y=1.01
)
plt.tight_layout()

out_fig = os.path.join(OUT_DIR, 'nihss_fc_within_cluster_3timepoints.png')
fig.savefig(out_fig, dpi=200, bbox_inches='tight')
plt.close()
print(f'\nSaved figure: {out_fig}')

# Save results CSV
res_df = pd.DataFrame(all_results)[['timepoint','network','cluster','n','rho','p','sig']]
csv_out = os.path.join(OUT_DIR, 'nihss_fc_within_cluster_3timepoints.csv')
res_df.to_csv(csv_out, index=False)
print(f'Saved CSV: {csv_out}')

print('\n── Summary ─────────────────────────────────────────────────────────')
print(res_df.to_string(index=False))
