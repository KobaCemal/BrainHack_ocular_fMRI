"""
1. Full parcel-level: for each disconnectome tract, Spearman r with each of
   400 Schaefer parcels' beta_r across subjects → spatial map per tract.
   FDR-BH across all tests.
2. Scatter plots: top 6 (one per tract, best parcel) + betar_Vis summaries.
"""
import sys
sys.path.insert(0, '/mnt/compneuro/comparison_of_methods/OMvar_paper')
import config
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
from nilearn import datasets as nlds
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import glob

INT_X      = Path(config.INT_DIR) / 'x_distance_analysis'
LESION_SRC = INT_X / 'lesion_masks' / 'NiftiLesions_original'
COPE_DIR   = Path(config.INT_DIR) / '2_Correlates' / 'allsessions_glm_zscored' / 'acute'
FIG_DIR    = Path(config.FIG_DIR) / 'x_distance_analysis'
FIG_DIR.mkdir(parents=True, exist_ok=True)
DISC_CSV   = '/mnt/compneuro/deepmreye_finetuning/disconnectome_tract_values.csv'

# ── 1. Reconstruct common_ids ─────────────────────────────────────────────────
lesion_files = sorted(glob.glob(str(LESION_SRC / 'sub-PAT*_A_lesion_ifhfix.nii.gz')))
lesion_ids   = [Path(f).name.replace('.nii.gz','').split('_A_')[0].replace('sub-','')
                for f in lesion_files]
cope_files   = sorted(glob.glob(str(COPE_DIR / 'PAT*_cope_fe.nii.gz')))
cope_set     = {Path(f).name.replace('.nii.gz','').replace('_cope_fe','') for f in cope_files}
common_ids   = [pid for pid in lesion_ids if pid in cope_set]
N            = len(common_ids)

# ── 2. Load Schaefer cope (N, 400) ────────────────────────────────────────────
cope_schaef = np.load(str(INT_X / 'cope_schaefer400.npy'))   # (N, 400)
atlas_obj     = nlds.fetch_atlas_schaefer_2018(n_rois=400, yeo_networks=7, resolution_mm=2)
parcel_labels = [l.decode() if isinstance(l, bytes) else l for l in atlas_obj.labels]
if parcel_labels[0].lower() in ('background', '0', ''):
    parcel_labels = parcel_labels[1:]
parcel_labels = parcel_labels[:400]
networks = [l.split('_')[2] for l in parcel_labels]

# ── 3. Load disconnectome ─────────────────────────────────────────────────────
disc = pd.read_csv(DISC_CSV)
disc['ID'] = disc['ID'].str.replace('sub-', '', regex=False)

tract_cols = ['disc_SLF_ipsi', 'disc_IFOF_ipsi', 'disc_ILF_ipsi']

betar_df = pd.DataFrame({'ID': common_ids})
merged   = betar_df.merge(disc[['ID'] + tract_cols], on='ID', how='inner')
n_merged = len(merged)
print(f'N merged: {n_merged}')

# reorder cope to match merged IDs
id_to_idx = {pid: i for i, pid in enumerate(common_ids)}
row_idx   = [id_to_idx[pid] for pid in merged['ID']]
cope_sub  = cope_schaef[row_idx, :]   # (n_merged, 400)

# ── 4. Restrict parcels to nb3 significant clusters ──────────────────────────
import nibabel as nib
from nilearn import image as nlimage
import sys
sys.path.insert(0, '/mnt/compneuro/comparison_of_methods/OMvar_paper')
INT_3 = Path(config.INT_DIR) / '3_TemporalFeatures'
neg_mask_img = nib.load(str(INT_3 / '3_neg_cluster_mask.nii.gz'))
pos_mask_img = nib.load(str(INT_3 / '3_pos_cluster_mask.nii.gz'))
atlas_obj2   = nlds.fetch_atlas_schaefer_2018(n_rois=400, yeo_networks=7, resolution_mm=2)
atlas_cl     = nlimage.resample_to_img(atlas_obj2.maps, neg_mask_img, interpolation='nearest')
atlas_cl_d   = atlas_cl.get_fdata().ravel()

def cluster_parcel_indices(mask_img, atlas_data):
    mask = mask_img.get_fdata().ravel() > 0
    pids = np.unique(atlas_data[mask & (atlas_data > 0)]).astype(int)
    return pids - 1   # 0-indexed

neg_idx = cluster_parcel_indices(neg_mask_img, atlas_cl_d)
pos_idx = cluster_parcel_indices(pos_mask_img, atlas_cl_d)
both_idx = np.unique(np.concatenate([neg_idx, pos_idx]))
print(f'Neg cluster parcels: {len(neg_idx)}, Pos cluster parcels: {len(pos_idx)}, Combined: {len(both_idx)}')

# ── 5. Full parcel-level correlations (cluster-restricted) ────────────────────
print('Running parcel-level correlations...')
rows = []
def zscore_outlier_mask(v, thresh=2.5):
    m, s = np.nanmean(v), np.nanstd(v)
    return np.abs(v - m) <= thresh * s if s > 0 else np.ones(len(v), dtype=bool)

for tc in tract_cols:
    x = merged[tc].values.astype(float)
    nonzero = x > 0
    print(f'{tc}: {nonzero.sum()} non-zero subjects')
    xz = (x - np.nanmean(x[nonzero])) / np.nanstd(x[nonzero])
    for j in both_idx:
        y = cope_sub[:, j]
        yz = (y - np.nanmean(y)) / np.nanstd(y)
        valid = (nonzero & np.isfinite(xz) & np.isfinite(yz)
                 & (np.abs(xz) <= 2.5) & (np.abs(yz) <= 2.5))
        if valid.sum() >= 10 and xz[valid].std() > 0:
            r, p = spearmanr(xz[valid], yz[valid])
            cluster = ('pos' if j in pos_idx else '') + ('neg' if j in neg_idx else '')
            rows.append({'tract': tc, 'parcel_idx': j,
                         'parcel_label': parcel_labels[j],
                         'network': networks[j],
                         'cluster': cluster,
                         'rho': r, 'pval': p, 'n': valid.sum()})

results = pd.DataFrame(rows)
print(f'Total tests: {len(results)}')

fin = results['pval'].notna()
results['qval'] = np.nan
_, qvals_fin, _, _ = multipletests(results.loc[fin, 'pval'], method='fdr_bh')
results.loc[fin, 'qval'] = qvals_fin

n_fdr = (results['qval'] < 0.05).sum()
n_unc = (results['pval'] < 0.05).sum()
print(f'FDR q<0.05: {n_fdr}  |  uncorrected p<0.05: {n_unc}')

out_path = str(INT_X / 'disconnectome_betar_parcel_results.csv')
results.to_csv(out_path, index=False)
print(f'Saved: {out_path}')

print('\nTop 20 by p-value:')
top20 = results.sort_values('pval').head(20)
print(top20[['tract', 'parcel_label', 'rho', 'pval', 'qval']].to_string(index=False))

# ── 5. Scatter plots: top parcel per tract ────────────────────────────────────
top_per_tract = results.sort_values('pval').groupby('tract').first().reset_index()
top_per_tract = top_per_tract.sort_values('pval')

fig, axes = plt.subplots(1, 3, figsize=(13, 4))
axes = axes.ravel()

for ax, (_, row) in zip(axes, top_per_tract.iterrows()):
    tc  = row['tract']
    j   = int(row['parcel_idx'])
    x   = merged[tc].values.astype(float)
    y   = cope_sub[:, j]
    xz  = (x - np.nanmean(x[x > 0])) / np.nanstd(x[x > 0])
    yz  = (y - np.nanmean(y)) / np.nanstd(y)
    valid = ((x > 0) & np.isfinite(xz) & np.isfinite(yz)
             & (np.abs(xz) <= 2.5) & (np.abs(yz) <= 2.5))
    ax.scatter(xz[valid], yz[valid], s=25, alpha=0.6, color='steelblue', edgecolors='none')
    m, b = np.polyfit(xz[valid], yz[valid], 1)
    xl = np.linspace(xz[valid].min(), xz[valid].max(), 100)
    ax.plot(xl, m * xl + b, color='crimson', lw=1.5)
    ax.set_xlabel(tc.replace('disc_', '').replace('_', ' ') + ' (z)', fontsize=9)
    ax.set_ylabel('beta_r (z)', fontsize=9)
    lbl = row['parcel_label'].replace('7Networks_', '')
    q_str = f"q={row['qval']:.3f}" if np.isfinite(row['qval']) else ''
    ax.set_title(f"{lbl}\nρ={row['rho']:+.3f}  p={row['pval']:.4f}  {q_str}", fontsize=8)

plt.suptitle('Top parcel per disconnectome tract × beta_r  (non-zero disconnection only)', fontsize=11)
plt.tight_layout()
fig_path = str(FIG_DIR / 'disconnectome_betar_scatter.png')
plt.savefig(fig_path, dpi=150)
plt.close()
print(f'\nSaved figure: {fig_path}')
