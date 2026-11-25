"""
Mass-univariate: for each Schaefer-400 parcel, Spearman correlation between
parcel-level lesion damage and parcel-level beta_r (oculomotor coupling) across
N=111 acute patients.
Lesion damage = mean binary lesion value in parcel voxels (fraction lesioned).
Beta_r = mean GLM cope value within parcel (cope_schaefer400.npy cached).
FDR-BH across 400 parcels.
"""
import sys
sys.path.insert(0, '/mnt/compneuro/comparison_of_methods/OMvar_paper')
import config
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
import glob
import nibabel as nib
from nilearn import datasets as nlds, image

INT_X      = Path(config.INT_DIR) / 'x_distance_analysis'
LESION_SRC = INT_X / 'lesion_masks' / 'NiftiLesions_original'
COPE_DIR   = Path(config.INT_DIR) / '2_Correlates' / 'allsessions_glm_zscored' / 'acute'

# ── 1. Reconstruct common_ids ─────────────────────────────────────────────────
lesion_files = sorted(glob.glob(str(LESION_SRC / 'sub-PAT*_A_lesion_ifhfix.nii.gz')))
lesion_ids   = [Path(f).name.replace('.nii.gz','').split('_A_')[0].replace('sub-','')
                for f in lesion_files]
cope_files   = sorted(glob.glob(str(COPE_DIR / 'PAT*_cope_fe.nii.gz')))
cope_set     = {Path(f).name.replace('.nii.gz','').replace('_cope_fe','') for f in cope_files}
common_ids   = [pid for pid in lesion_ids if pid in cope_set]
N            = len(common_ids)
print(f'N = {N}')

# ── 2. Load lesion_flat and build parcel-level lesion damage ──────────────────
d           = np.load(str(INT_X / 'lesion_array_4mm.npz'), allow_pickle=True)
lesion_flat = d['lesion_flat']   # (N, n_vox) full volume raveled
lesion_bin  = (lesion_flat > 0.5).astype(np.float32)

ref_img     = nib.load(str(COPE_DIR / f'{common_ids[0]}_cope_fe.nii.gz'))
ref_shape   = ref_img.shape[:3]
print(f'Reference space: {ref_shape}  n_vox: {lesion_flat.shape[1]}')

atlas_obj  = nlds.fetch_atlas_schaefer_2018(n_rois=400, yeo_networks=7, resolution_mm=2)
atlas_4mm  = image.resample_to_img(atlas_obj.maps, ref_img, interpolation='nearest')
atlas_data = atlas_4mm.get_fdata(dtype=np.float32).ravel()   # (n_vox,)

parcel_ids = np.arange(1, 401)
lesion_parc = np.zeros((N, 400), dtype=np.float32)   # fraction lesioned per parcel
for j, pid_roi in enumerate(parcel_ids):
    mask = (atlas_data == pid_roi)
    if mask.sum() > 0:
        lesion_parc[:, j] = lesion_bin[:, mask].mean(axis=1)
    else:
        lesion_parc[:, j] = np.nan

n_empty = np.isnan(lesion_parc[0]).sum()
print(f'Empty parcels: {n_empty}/400')
print(f'Lesion damage range: {lesion_parc[~np.isnan(lesion_parc)].min():.4f} – {lesion_parc[~np.isnan(lesion_parc)].max():.4f}')

# ── 3. Load cached Schaefer beta_r ───────────────────────────────────────────
cope_schaef = np.load(str(INT_X / 'cope_schaefer400.npy'))   # (N, 400)
print(f'cope_schaefer shape: {cope_schaef.shape}')

# ── 4. Spearman correlation per parcel ────────────────────────────────────────
rhos = np.full(400, np.nan)
pvals = np.full(400, np.nan)

for j in range(400):
    x = lesion_parc[:, j]
    y = cope_schaef[:, j]
    valid = np.isfinite(x) & np.isfinite(y)
    if valid.sum() >= 10:
        r, p = spearmanr(x[valid], y[valid])
        rhos[j]  = r
        pvals[j] = p

n_tested = np.isfinite(rhos).sum()
print(f'Parcels tested: {n_tested}/400')

# ── 5. FDR-BH ─────────────────────────────────────────────────────────────────
valid_idx = np.where(np.isfinite(pvals))[0]
_, qvals_valid, _, _ = multipletests(pvals[valid_idx], method='fdr_bh')
qvals = np.full(400, np.nan)
qvals[valid_idx] = qvals_valid

n_sig = (qvals < 0.05).sum()
print(f'FDR q<0.05: {n_sig} parcels')
print(f'Uncorrected p<0.05: {(pvals < 0.05).sum()} parcels')

# ── 6. Get parcel labels ───────────────────────────────────────────────────────
parcel_labels = [l.decode() if isinstance(l, bytes) else l for l in atlas_obj.labels]
# skip background if present
if parcel_labels[0].lower() in ('background', '0', ''):
    parcel_labels = parcel_labels[1:]

# ── 7. Save results ───────────────────────────────────────────────────────────
results = pd.DataFrame({
    'parcel_id':    parcel_ids,
    'parcel_label': parcel_labels[:400],
    'rho':          rhos,
    'pval':         pvals,
    'qval':         qvals,
})

out_path = str(INT_X / 'massunivariate_lesion_betar.csv')
results.to_csv(out_path, index=False)
print(f'Saved: {out_path}')

# ── 8. Top results ────────────────────────────────────────────────────────────
top = results.dropna().sort_values('pval').head(20)
print('\nTop 20 by p-value:')
print(top[['parcel_label', 'rho', 'pval', 'qval']].to_string(index=False))
