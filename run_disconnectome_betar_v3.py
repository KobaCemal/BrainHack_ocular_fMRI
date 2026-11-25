"""
Correlate disconnectome tract values (SLF/IFOF/ILF ipsi) with mean beta_r
within the two nb3 clusters (positive=visual, negative=prefrontal).
2 outcomes × 3 tracts = 6 tests. z-score + |z|>2.5 outlier removal.
"""
import sys
sys.path.insert(0, '/mnt/compneuro/comparison_of_methods/OMvar_paper')
import config
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
import nibabel as nib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import glob

INT_X      = Path(config.INT_DIR) / 'x_distance_analysis'
INT_3      = Path(config.INT_DIR) / '3_TemporalFeatures'
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
print(f'N = {N}')

# ── 2. Load cluster masks ─────────────────────────────────────────────────────
pos_img  = nib.load(str(INT_3 / '3_pos_cluster_mask.nii.gz'))
neg_img  = nib.load(str(INT_3 / '3_neg_cluster_mask.nii.gz'))
pos_mask = pos_img.get_fdata().ravel() > 0
neg_mask = neg_img.get_fdata().ravel() > 0
print(f'Pos cluster: {pos_mask.sum()} voxels  |  Neg cluster: {neg_mask.sum()} voxels')

# ── 3. Compute mean cope within each cluster per subject ──────────────────────
betar_pos = np.full(N, np.nan)
betar_neg = np.full(N, np.nan)

for i, pid in enumerate(common_ids):
    cope_path = str(COPE_DIR / f'{pid}_cope_fe.nii.gz')
    img = nib.load(cope_path)
    d   = img.get_fdata(dtype=np.float32).ravel()
    betar_pos[i] = d[pos_mask].mean()
    betar_neg[i] = d[neg_mask].mean()

print(f'betar_pos: mean={np.nanmean(betar_pos):.4f}  std={np.nanstd(betar_pos):.4f}')
print(f'betar_neg: mean={np.nanmean(betar_neg):.4f}  std={np.nanstd(betar_neg):.4f}')

# ── 4. Load disconnectome ─────────────────────────────────────────────────────
disc = pd.read_csv(DISC_CSV)
disc['ID'] = disc['ID'].str.replace('sub-', '', regex=False)
tract_cols = ['disc_SLF_ipsi', 'disc_IFOF_ipsi', 'disc_ILF_ipsi']

betar_df = pd.DataFrame({'ID': common_ids, 'betar_pos': betar_pos, 'betar_neg': betar_neg})
merged   = betar_df.merge(disc[['ID'] + tract_cols], on='ID', how='inner')
print(f'Merged: {len(merged)}')

# ── 5. z-score + outlier removal, Spearman ────────────────────────────────────
def zscored(v):
    m, s = np.nanmean(v), np.nanstd(v)
    return (v - m) / s if s > 0 else np.zeros_like(v)

rows = []
for tc in tract_cols:
    x     = merged[tc].values.astype(float)
    nonzero = x > 0
    xz    = zscored(x[nonzero])  # z-score within non-zero subjects
    xz_full = np.full(len(x), np.nan)
    xz_full[nonzero] = xz

    for outcome, label in [('betar_pos', 'pos_cluster'), ('betar_neg', 'neg_cluster')]:
        y  = merged[outcome].values.astype(float)
        yz = zscored(y)
        valid = (nonzero & np.isfinite(xz_full) & np.isfinite(yz)
                 & (np.abs(xz_full) <= 2.5) & (np.abs(yz) <= 2.5))
        n = valid.sum()
        if n >= 10:
            r, p = spearmanr(xz_full[valid], yz[valid])
            rows.append({'tract': tc, 'outcome': label,
                         'rho': r, 'pval': p, 'n': n,
                         'xz': xz_full, 'yz': yz, 'valid': valid})

results = pd.DataFrame([{k: v for k, v in r.items() if k not in ('xz','yz','valid')}
                        for r in rows])
_, qvals, _, _ = multipletests(results['pval'], method='fdr_bh')
results['qval'] = qvals

print(f'\nResults ({len(results)} tests):')
print(results[['tract','outcome','rho','pval','qval','n']].to_string(index=False))

out_path = str(INT_X / 'disconnectome_betar_cluster_results.csv')
results.to_csv(out_path, index=False)
print(f'\nSaved: {out_path}')

# ── 6. Scatter plots ──────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 3, figsize=(13, 8))
for col_i, tc in enumerate(tract_cols):
    for row_i, (outcome, label) in enumerate([('betar_pos', 'Pos cluster (visual)'),
                                               ('betar_neg', 'Neg cluster (prefrontal)')]):
        ax  = axes[row_i, col_i]
        row = next(r for r in rows if r['tract'] == tc and r['outcome'] == label.split(' ')[0].lower() + '_cluster')
        xv  = row['xz'][row['valid']]
        yv  = row['yz'][row['valid']]
        ax.scatter(xv, yv, s=25, alpha=0.6, color='steelblue', edgecolors='none')
        m, b = np.polyfit(xv, yv, 1)
        xl = np.linspace(xv.min(), xv.max(), 100)
        ax.plot(xl, m * xl + b, color='crimson', lw=1.5)
        res_row = results[(results['tract'] == tc) & (results['outcome'] == row['outcome'])].iloc[0]
        ax.set_title(f"{label}\nρ={res_row['rho']:+.3f}  p={res_row['pval']:.4f}  q={res_row['qval']:.3f}", fontsize=9)
        ax.set_xlabel(tc.replace('disc_','').replace('_',' ') + ' (z)', fontsize=9)
        ax.set_ylabel('beta_r (z)', fontsize=9)

plt.suptitle('Disconnectome ipsi × nb3 cluster beta_r  (non-zero, |z|≤2.5)', fontsize=11)
plt.tight_layout()
fig_path = str(FIG_DIR / 'disconnectome_cluster_betar_scatter.png')
plt.savefig(fig_path, dpi=150)
plt.close()
print(f'Saved figure: {fig_path}')
