"""
NIHSS cluster × Network FC boxplots — 2 weeks, 3 months, 1 year
-----------------------------------------------------------------
For each of the 7 Yeo networks and each of 3 timepoints, compute
mean within-network FC from mean_corrmats and group by session-
specific nih_total into 4 NIHSS severity clusters.

One combined figure: 7 panels (one per network), grouped boxplots
with NIHSS cluster on x-axis and timepoint as hue.
"""

import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
import os, warnings
warnings.filterwarnings('ignore')

# ── Paths ───────────────────────────────────────────────────────────────────
DATA_PATH = '/home/cemal/Desktop/Opus/stroke_clinical_and_neural_data_with_lesion_disconnectome.pkl'
MAP_PATH  = '/home/cemal/Desktop/Opus/R4_analyses/schaefer400_network_mapping.csv'
OUT_DIR   = '/home/cemal/Desktop/Opus/NIHSS_FC_correlation'

# ── Load ────────────────────────────────────────────────────────────────────
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
net_to_idx = {net: net_map[net_map['network'] == net]['roi_idx'].values for net in networks}

# ── Filter stroke patients ───────────────────────────────────────────────────
stroke = data[data['subj_type'] == 0].copy()
excluded = stroke.loc[stroke['inclusion_notes'].notna(), 'ID'].unique()
stroke = stroke[~stroke['ID'].isin(excluded)].copy()

# ── Session mapping ──────────────────────────────────────────────────────────
sessions = [
    ('acute',     '2 weeks'),
    ('followup',  '3 months'),
    ('followup2', '1 year'),
]

# ── NIHSS cluster settings ───────────────────────────────────────────────────
cat_bins   = [-0.1, 1, 5, 10, 100]
cat_labels = ['0–1', '2–5', '5–10', '10+']

# ── Helper: within-network mean FC ──────────────────────────────────────────
def within_network_mean(mat, idx):
    sub = mat[np.ix_(idx, idx)]
    return float(np.nanmean(sub[np.triu_indices(len(idx), k=1)]))

# ── Build combined long-format DataFrame ────────────────────────────────────
all_rows = []
for sess_key, sess_label in sessions:
    sess_df = stroke[stroke['Session'] == sess_key].copy()
    for _, row in sess_df.iterrows():
        mat   = row['mean_corrmats']
        nihss = row['nih_total']
        if not isinstance(mat, np.ndarray) or mat.shape != (400, 400):
            continue
        if pd.isna(nihss):
            continue
        rec = {
            'ID':        row['ID'],
            'timepoint': sess_label,
            'nih_total': nihss,
            'nihss_cat': pd.cut([nihss], bins=cat_bins, labels=cat_labels)[0],
        }
        for net in networks:
            rec[f'fc_{net}'] = within_network_mean(mat, net_to_idx[net])
        all_rows.append(rec)

df = pd.DataFrame(all_rows)
df['nihss_cat'] = pd.Categorical(df['nihss_cat'], categories=cat_labels, ordered=True)
df['timepoint'] = pd.Categorical(df['timepoint'],
                                  categories=['2 weeks', '3 months', '1 year'], ordered=True)

print('Sample sizes per timepoint:')
print(df.groupby('timepoint', observed=True)['ID'].count().to_string())
print('\nNIHSS cluster × timepoint counts:')
print(df.groupby(['timepoint','nihss_cat'], observed=True)['ID'].count().unstack().to_string())

# ── Plot settings ────────────────────────────────────────────────────────────
tp_colors  = {'2 weeks': '#3498db', '3 months': '#e67e22', '1 year': '#2ecc71'}
tp_order   = ['2 weeks', '3 months', '1 year']
n_clusters = len(cat_labels)
n_tp       = len(tp_order)

# Positions: clusters spaced 1 unit apart, timepoints offset within each cluster
cluster_gap  = 1.0
tp_width     = 0.22
tp_offsets   = np.linspace(-(n_tp - 1) * tp_width / 2,
                            (n_tp - 1) * tp_width / 2,
                            n_tp)
cluster_centers = np.arange(n_clusters) * cluster_gap

# ── Main figure: 7 panels ────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 4, figsize=(26, 12))
axes = axes.flatten()

for ax_i, net in enumerate(networks):
    ax   = axes[ax_i]
    fc_c = f'fc_{net}'

    for ci, cat in enumerate(cat_labels):
        for ti, tp in enumerate(tp_order):
            vals = df[(df['nihss_cat'] == cat) & (df['timepoint'] == tp)][fc_c].dropna().values
            if len(vals) == 0:
                continue

            xpos = cluster_centers[ci] + tp_offsets[ti]
            bp = ax.boxplot(
                vals,
                positions=[xpos],
                widths=tp_width * 0.85,
                patch_artist=True,
                notch=False,
                showfliers=True,
                flierprops=dict(marker='o', markersize=3,
                                markerfacecolor=tp_colors[tp], alpha=0.5,
                                markeredgewidth=0.3),
                medianprops=dict(color='black', linewidth=1.5),
                boxprops=dict(facecolor=tp_colors[tp], alpha=0.75,
                              linewidth=0.8),
                whiskerprops=dict(linewidth=0.8),
                capprops=dict(linewidth=0.8),
            )

            # n label below each box
            ax.text(xpos, ax.get_ylim()[0] if ax_i == 0 else -999,
                    f'n={len(vals)}', ha='center', va='top', fontsize=5.5,
                    color='#555')

    # x-axis ticks at cluster centres
    ax.set_xticks(cluster_centers)
    ax.set_xticklabels(cat_labels, fontsize=9)
    ax.set_xlim(cluster_centers[0] - 0.55, cluster_centers[-1] + 0.55)
    ax.set_title(net_labels[net], fontsize=11, fontweight='bold', pad=5)
    ax.set_xlabel('NIHSS cluster', fontsize=9)
    ax.set_ylabel('Mean within-network FC (r)', fontsize=9)
    ax.tick_params(labelsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Cluster separator lines
    for c in cluster_centers[:-1]:
        ax.axvline(c + 0.5, color='#ccc', lw=0.8, ls=':')

# Fix the n= labels (replot them after ylim is stable)
for ax_i, net in enumerate(networks):
    ax   = axes[ax_i]
    fc_c = f'fc_{net}'
    ymin = ax.get_ylim()[0]
    for ci, cat in enumerate(cat_labels):
        for ti, tp in enumerate(tp_order):
            vals = df[(df['nihss_cat'] == cat) & (df['timepoint'] == tp)][fc_c].dropna().values
            if len(vals) == 0:
                continue
            xpos = cluster_centers[ci] + tp_offsets[ti]
            ax.text(xpos, ymin, f'n={len(vals)}',
                    ha='center', va='top', fontsize=5.5, color='#555')

# Remove unused panel
axes[-1].set_visible(False)

# Legend
patches = [mpatches.Patch(color=tp_colors[tp], label=tp, alpha=0.8) for tp in tp_order]
fig.legend(handles=patches, title='Timepoint', loc='lower right',
           bbox_to_anchor=(0.97, 0.05), fontsize=11, title_fontsize=12,
           frameon=True, edgecolor='grey')

fig.suptitle(
    'Mean within-network FC by NIHSS severity cluster\nat 2 weeks, 3 months, and 1 year post-stroke',
    fontsize=13, fontweight='bold', y=1.01
)

plt.tight_layout()
out_path = os.path.join(OUT_DIR, 'nihss_fc_boxplots_3timepoints.png')
fig.savefig(out_path, dpi=200, bbox_inches='tight')
plt.close()
print(f'\nSaved: {out_path}')

# ── Save long-format data table ──────────────────────────────────────────────
fc_cols = [f'fc_{n}' for n in networks]
out_df  = df[['ID','timepoint','nih_total','nihss_cat'] + fc_cols].copy()
out_df.columns = ['ID','timepoint','nih_total','nihss_cat'] + [net_labels[n] for n in networks]
csv_path = os.path.join(OUT_DIR, 'nihss_fc_boxplot_data.csv')
out_df.to_csv(csv_path, index=False)
print(f'Saved: {csv_path}')
