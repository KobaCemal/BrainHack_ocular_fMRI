"""
Correlate disconnectome tract values (SLF/IFOF/ILF ipsi/contra) with
mean beta_r (oculomotor coupling cope) whole-brain and per 7-network,
across acute patients with both datasets.
Spearman r, FDR-BH across all tests.
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
import glob

INT_X      = Path(config.INT_DIR) / 'x_distance_analysis'
LESION_SRC = INT_X / 'lesion_masks' / 'NiftiLesions_original'
COPE_DIR   = Path(config.INT_DIR) / '2_Correlates' / 'allsessions_glm_zscored' / 'acute'
DISC_CSV   = '/mnt/compneuro/deepmreye_finetuning/disconnectome_tract_values.csv'

# ── 1. Reconstruct common_ids ─────────────────────────────────────────────────
lesion_files = sorted(glob.glob(str(LESION_SRC / 'sub-PAT*_A_lesion_ifhfix.nii.gz')))
lesion_ids   = [Path(f).name.replace('.nii.gz','').split('_A_')[0].replace('sub-','')
                for f in lesion_files]
cope_files   = sorted(glob.glob(str(COPE_DIR / 'PAT*_cope_fe.nii.gz')))
cope_set     = {Path(f).name.replace('.nii.gz','').replace('_cope_fe','') for f in cope_files}
common_ids   = [pid for pid in lesion_ids if pid in cope_set]
N            = len(common_ids)
print(f'beta_r subjects: {N}')

# ── 2. Load Schaefer cope and compute mean per subject per network ─────────────
cope_schaef = np.load(str(INT_X / 'cope_schaefer400.npy'))   # (N, 400)

atlas_obj     = nlds.fetch_atlas_schaefer_2018(n_rois=400, yeo_networks=7, resolution_mm=2)
parcel_labels = [l.decode() if isinstance(l, bytes) else l for l in atlas_obj.labels]
if parcel_labels[0].lower() in ('background', '0', ''):
    parcel_labels = parcel_labels[1:]
parcel_labels = parcel_labels[:400]

networks      = [l.split('_')[2] for l in parcel_labels]
network_names = sorted(set(networks))
print(f'Networks: {network_names}')

# whole-brain mean and per-network mean
betar_df = pd.DataFrame({'ID': common_ids})
betar_df['betar_whole'] = np.nanmean(cope_schaef, axis=1)
for net in network_names:
    idx = [i for i, n in enumerate(networks) if n == net]
    betar_df[f'betar_{net}'] = np.nanmean(cope_schaef[:, idx], axis=1)

print(f'beta_r columns: {[c for c in betar_df.columns if c != "ID"]}')

# ── 3. Load disconnectome ─────────────────────────────────────────────────────
disc = pd.read_csv(DISC_CSV)
disc['ID'] = disc['ID'].str.replace('sub-', '', regex=False)
print(f'Disconnectome subjects: {len(disc)}')

tract_cols = ['disc_SLF_ipsi', 'disc_SLF_contra',
              'disc_IFOF_ipsi', 'disc_IFOF_contra',
              'disc_ILF_ipsi',  'disc_ILF_contra']

# ── 4. Merge ──────────────────────────────────────────────────────────────────
merged = betar_df.merge(disc[['ID'] + tract_cols], on='ID', how='inner')
n_merged = len(merged)
print(f'Merged (both datasets): {n_merged}')

# ── 5. Spearman correlations ──────────────────────────────────────────────────
betar_cols = [c for c in betar_df.columns if c != 'ID']
rows = []
for tc in tract_cols:
    for bc in betar_cols:
        x = merged[tc].values
        y = merged[bc].values
        valid = np.isfinite(x) & np.isfinite(y)
        if valid.sum() >= 10:
            r, p = spearmanr(x[valid], y[valid])
            rows.append({'tract': tc, 'betar': bc, 'rho': r, 'pval': p, 'n': valid.sum()})

results = pd.DataFrame(rows)

# FDR-BH (only on finite p-values)
results['qval'] = np.nan
fin = results['pval'].notna()
if fin.sum() > 0:
    _, qvals_fin, _, _ = multipletests(results.loc[fin, 'pval'], method='fdr_bh')
    results.loc[fin, 'qval'] = qvals_fin

n_sig_fdr = (results['qval'] < 0.05).sum()
n_sig_unc = (results['pval'] < 0.05).sum()
print(f'\nTotal tests: {len(results)}')
print(f'FDR q<0.05: {n_sig_fdr}')
print(f'Uncorrected p<0.05: {n_sig_unc}')

# ── 6. Report ─────────────────────────────────────────────────────────────────
top = results.sort_values('pval').head(20)
print('\nTop 20:')
print(top[['tract', 'betar', 'rho', 'pval', 'qval', 'n']].to_string(index=False))

out_path = str(INT_X / 'disconnectome_betar_results.csv')
results.to_csv(out_path, index=False)
print(f'\nSaved: {out_path}')
