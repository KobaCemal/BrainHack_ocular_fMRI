"""Run this script once to create scratch_1_Eye_Brain_lag.ipynb."""
import json, os

def code(src):
    return {"cell_type": "code", "execution_count": None,
            "metadata": {}, "outputs": [], "source": src}

def md(src):
    return {"cell_type": "markdown", "metadata": {}, "source": src}

cells = []

# ── Cell 0: title ─────────────────────────────────────────────────────────────
cells.append(md(
"""# Notebook 1: Temporal structure of eye–BOLD coupling

**Question:** What is the temporal structure of the coupling between spontaneous \
eye movements (r_trtrcorr) and BOLD signal across the visual and oculomotor hierarchy?

**Approach:** Pearson cross-correlation between Nudelman-residualized ROI-mean BOLD \
and r_trtrcorr, at lags −5 to +5 TRs (−10 s to +10 s).

**Lag convention:** positive lag = eye movement *precedes* BOLD. \
A peak at lag +2 (+4 s) means the eye signal 4 s ago best predicts current BOLD.

**Pipeline:** Nudelman no-BP — 16 confounds, AR1, 6 mm FWHM, no bandpass filter. \
Defined in `config.py`."""
))

# ── Cell 1: imports ───────────────────────────────────────────────────────────
cells.append(code(
"""import os, sys, pickle
sys.path.insert(0, '/mnt/compneuro/comparison_of_methods/OMvar_paper')
from config import *

import numpy as np
import nibabel as nib
import pandas as pd
from scipy.stats import pearsonr, sem as scipy_sem
from scipy.ndimage import gaussian_filter
from joblib import Parallel, delayed
%matplotlib inline
import matplotlib.pyplot as plt
from nilearn import datasets, image

os.makedirs(os.path.join(FIG_DIR, '1_Eye_Brain_lag'), exist_ok=True)
os.makedirs(INT_DIR, exist_ok=True)

meta = pd.read_pickle(PKL)
groups, q25, q75 = get_groups(meta)
controls = groups['Controls']

print(f"PKL: {meta.shape}")
for name, df in groups.items():
    print(f"  {name}: {df['ID'].nunique()} subjects ({len(df)} rows)")
print(f"  Lesion quartile cutoffs — Q25: {q25:.0f} mm³  Q75: {q75:.0f} mm³")"""
))

# ── Cell 2: ROI definitions (markdown) ────────────────────────────────────────
cells.append(md(
"""## ROI definitions

| ROI | Source | Label(s) / coordinates |
|-----|--------|------------------------|
| LGN  | Juelich | GM Lateral geniculate body |
| V1   | Juelich | GM Visual cortex V1 BA17 |
| V2   | Juelich | GM Visual cortex V2 BA18 |
| V3V  | Juelich | GM Visual cortex V3V |
| V4   | Juelich | GM Visual cortex V4 |
| V5   | Juelich | GM Visual cortex V5 |
| IPS  | Juelich | GM Anterior intra-parietal sulcus hIP1/2/3 |
| FEF  | Sphere  | bilateral ±28,−6,+50 mm  r=8 mm (Amiez & Petrides 2009) |
| SEF  | Sphere  | 0,+4,+52 mm  r=8 mm |

Atlas: Juelich maxprob-thr0-1mm, resampled to MNI152NLin2009cAsym 2 mm \
(fMRIPrep output space). Masks saved to `intermediates/` on first run."""
))

# ── Cell 3: build ROI masks ───────────────────────────────────────────────────
cells.append(code(
"""MASK_CACHE    = os.path.join(INT_DIR, '1_roi_masks.npz')
AFFINE_CACHE  = os.path.join(INT_DIR, '1_ref_affine.npy')
SHAPE_CACHE   = os.path.join(INT_DIR, '1_ref_shape.npy')

if os.path.exists(MASK_CACHE):
    d          = np.load(MASK_CACHE)
    roi_masks  = {k: d[k].astype(bool) for k in d.files}
    ref_affine = np.load(AFFINE_CACHE)
    ref_shape  = tuple(np.load(SHAPE_CACHE).astype(int))
    print("ROI masks loaded from cache.")
else:
    # Pick reference image from first available control run
    first_row = controls.iloc[0]
    first_sub = first_row['ID']
    first_run = int(str(first_row.get('Run', '1')).split('-')[0])
    ref_path  = os.path.join(
        FMRIPREP, f'sub-{first_sub}', 'ses-control', 'func',
        f'sub-{first_sub}_ses-control_task-rest_run-{first_run}'
        '_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz')
    ref_img    = nib.load(ref_path)
    ref_affine = ref_img.affine
    ref_shape  = ref_img.shape[:3]

    juelich      = datasets.fetch_atlas_juelich('maxprob-thr0-1mm')
    atlas_img    = image.load_img(juelich.maps)   # accepts path or image object
    atlas_res    = image.resample_img(
        atlas_img, target_affine=ref_affine,
        target_shape=ref_shape, interpolation='nearest')
    atlas_data   = atlas_res.get_fdata().astype(int)
    atlas_labels = list(juelich.labels)

    roi_masks = {}
    for roi, label_list in JUELICH_LABELS.items():
        mask = np.zeros(ref_shape, dtype=bool)
        for lbl in label_list:
            hits = [i for i, l in enumerate(atlas_labels) if l == lbl]
            if not hits:
                print(f'  WARNING: label not found in atlas: "{lbl}"')
            for idx in hits:
                mask |= (atlas_data == idx)
        roi_masks[roi] = mask

    # Sphere ROIs: compute MNI coords for every voxel once
    xi, yi, zi  = np.mgrid[:ref_shape[0], :ref_shape[1], :ref_shape[2]]
    vox_homo    = np.vstack([xi.ravel(), yi.ravel(), zi.ravel(),
                             np.ones(xi.size)])
    mni_coords  = (ref_affine @ vox_homo)[:3].T   # (n_vox, 3)

    for roi, params in SPHERE_ROIS.items():
        mask = np.zeros(ref_shape, dtype=bool)
        for center in params['centers']:
            dist  = np.sqrt(((mni_coords - np.array(center))**2).sum(axis=1))
            mask |= (dist <= params['radius']).reshape(ref_shape)
        roi_masks[roi] = mask

    np.savez(MASK_CACHE, **roi_masks)
    np.save(AFFINE_CACHE, ref_affine)
    np.save(SHAPE_CACHE,  np.array(ref_shape))
    print("ROI masks built and cached.")

print("\\nROI voxel counts:")
for roi in ROI_ORDER:
    n = roi_masks.get(roi, np.array([])).sum()
    print(f"  {roi:5s}: {n:5d} voxels")"""
))

# ── Cell 4: glassbrain ───────────────────────────────────────────────────────
cells.append(code(
"""import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from nilearn import plotting

# Combined label image: each ROI gets a unique integer
combined = np.zeros(ref_shape, dtype=np.float32)
for i, roi in enumerate(ROI_ORDER, start=1):
    combined[roi_masks[roi]] = i
combined_img = nib.Nifti1Image(combined, ref_affine)

cmap_combined = mcolors.ListedColormap([ROI_COLORS[roi] for roi in ROI_ORDER])

display = plotting.plot_glass_brain(
    combined_img,
    display_mode='lzr',
    colorbar=False,
    annotate=True,
    cmap=cmap_combined,
    vmin=0.5, vmax=len(ROI_ORDER) + 0.5,
    alpha=0.9,
    title='ROI masks — MNI152NLin2009cAsym'
)

fig_gb = plt.gcf()
handles = [mpatches.Patch(color=ROI_COLORS[roi], label=roi) for roi in ROI_ORDER]
fig_gb.legend(handles=handles, loc='lower center', fontsize=9,
              ncol=len(ROI_ORDER), bbox_to_anchor=(0.5, -0.02))
fig_gb.set_size_inches(14, 4)

OUT_GB = os.path.join(FIG_DIR, '1_Eye_Brain_lag', 'roi_glassbrain')
fig_gb.savefig(OUT_GB + '.png', dpi=150, bbox_inches='tight')
print(f'Saved: {OUT_GB}.png')
plt.show()"""
))

# ── Cell 5: method (markdown) ─────────────────────────────────────────────────
cells.append(md(
r"""## Cross-correlation method

For each control subject, for each run:

1. Load fMRI, trim 1 TR each end, smooth 6 mm FWHM
2. Nudelman 16-param + linear trend confound regression via lstsq
3. Extract mean residual BOLD within each ROI mask
4. For each lag $\ell \in [-5, \ldots, +5]$ TRs:

$$r(\ell) = \operatorname{Pearson}\!\bigl(\mathrm{BOLD}(t),\; \mathrm{eye}(t - \ell)\bigr)$$

where $\mathrm{eye}(t-\ell)$ means the eye signal is shifted **forward** by $\ell$ samples
(eye leads BOLD by $\ell \times 2$ s when $\ell > 0$).

5. Correlation curves are averaged across runs per subject, then group mean ± SEM across subjects.

**Eye movement predictor:** `r_trtrcorr` (NaN → 0 for missing TRs)."""
))

# ── Cell 5: compute cross-correlations ───────────────────────────────────────
cells.append(code(
"""LAGS        = list(range(-10, 11))  # -20 s to +20 s
XCORR_CACHE = os.path.join(INT_DIR, '1_xcorr_controls.pkl')

def _subject_xcorr(row, roi_masks, lags, session='control'):
    sub_id  = row['ID']
    eye_arr = row['r_trtrcorr']
    if not isinstance(eye_arr, np.ndarray) or eye_arr.ndim < 2:
        return sub_id, None

    func_dir = os.path.join(FMRIPREP, f'sub-{sub_id}', f'ses-{session}', 'func')
    if not os.path.isdir(func_dir):
        return sub_id, None

    run_corrs = {roi: [] for roi in roi_masks}

    for col_idx, run_num in parse_runs(row):
        fp = os.path.join(func_dir,
            f'sub-{sub_id}_ses-{session}_task-rest_run-{run_num}'
            '_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz')
        cp = os.path.join(func_dir,
            f'sub-{sub_id}_ses-{session}_task-rest_run-{run_num}'
            '_desc-confounds_timeseries.tsv')
        if not os.path.exists(fp) or not os.path.exists(cp):
            continue
        if col_idx >= eye_arr.shape[1]:
            continue
        eye_ts = eye_arr[:, col_idx]
        if np.isnan(eye_ts).all() or np.isfinite(eye_ts).sum() < 20:
            continue

        img  = nib.load(fp)
        data = img.get_fdata(dtype=np.float32)[:, :, :, 1:-1]
        n    = data.shape[3]

        vox_mm = np.sqrt((img.affine[:3, :3]**2).sum(axis=0))
        sigma  = FWHM_MM / (2.3548 * vox_mm)
        for t in range(n):
            data[:, :, :, t] = gaussian_filter(data[:, :, :, t], sigma=sigma)
        data2d = data.reshape(-1, n).T.astype(np.float64)

        conf_df = pd.read_csv(cp, sep='\\t').iloc[1:-1].reset_index(drop=True)
        C = np.zeros((n, len(NUDELMAN_COLS)), dtype=np.float64)
        for ci, col in enumerate(NUDELMAN_COLS):
            if col in conf_df.columns:
                C[:, ci] = np.nan_to_num(
                    conf_df[col].values.astype(float)[:n], nan=0.0)
        lin       = np.arange(n, dtype=np.float64); lin -= lin.mean()
        confounds = np.column_stack([np.ones(n), C, lin])
        betas, _, _, _ = np.linalg.lstsq(confounds, data2d, rcond=None)
        resid = data2d - confounds @ betas   # (n_trs, n_vox)

        eye_c = np.nan_to_num(eye_ts[:n], nan=0.0)
        eye_c -= eye_c.mean()
        if eye_c.std() < 1e-10:
            continue

        for roi, mask in roi_masks.items():
            mf = mask.flatten()
            if mf.sum() == 0:
                continue
            bold = resid[:, mf].mean(axis=1)
            bold_dm = bold - bold.mean()
            bold_sd = bold_dm.std()
            if bold_sd < 1e-10:
                continue

            lag_rs = []
            for lag in lags:
                shifted = np.zeros(n)
                if lag == 0:
                    shifted = eye_c.copy()
                elif lag > 0:
                    shifted[lag:] = eye_c[:n - lag]
                else:
                    shifted[:n + lag] = eye_c[-lag:]
                sh_sd = shifted.std()
                if sh_sd < 1e-10:
                    lag_rs.append(0.0)
                else:
                    r = np.dot(bold_dm, shifted) / (n * bold_sd * sh_sd)
                    lag_rs.append(float(np.clip(r, -1, 1)))
            run_corrs[roi].append(lag_rs)

    result = {
        roi: np.nanmean(runs, axis=0) if runs else np.full(len(lags), np.nan)
        for roi, runs in run_corrs.items()
    }
    return sub_id, result

if os.path.exists(XCORR_CACHE):
    with open(XCORR_CACHE, 'rb') as f:
        xcorr_all = pickle.load(f)
    print(f"Loaded from cache: {len(xcorr_all)} subjects")
else:
    rows_list = [row for _, row in controls.iterrows()]
    results   = Parallel(n_jobs=-1, prefer='processes', verbose=2)(
        delayed(_subject_xcorr)(row, roi_masks, LAGS) for row in rows_list
    )
    xcorr_all = {sid: res for sid, res in results if res is not None}
    with open(XCORR_CACHE, 'wb') as f:
        pickle.dump(xcorr_all, f)
    print(f"Computed and cached: {len(xcorr_all)} subjects")"""
))

# ── Cell 6: stats table ───────────────────────────────────────────────────────
cells.append(code(
"""lags_s = np.array(LAGS) * TR   # lags in seconds

subj_curves = {roi: [] for roi in ROI_ORDER}
for sub_id, res in xcorr_all.items():
    for roi in ROI_ORDER:
        if roi in res and not np.isnan(res[roi]).all():
            subj_curves[roi].append(res[roi])

means = {roi: np.nanmean(subj_curves[roi], axis=0)
         if subj_curves[roi] else np.full(len(LAGS), np.nan)
         for roi in ROI_ORDER}
sems  = {roi: scipy_sem(subj_curves[roi], axis=0, nan_policy='omit')
         if len(subj_curves[roi]) > 1 else np.full(len(LAGS), np.nan)
         for roi in ROI_ORDER}

lag0_idx = LAGS.index(0)

print(f"{'ROI':5s}  {'n':4s}  {'peak_lag_TR':11s}  {'peak_lag_s':10s}  "
      f"{'peak_r':8s}  {'r_lag0':8s}")
print('-' * 58)
peak_lag_TR = {}
for roi in ROI_ORDER:
    m = means[roi]
    n = len(subj_curves[roi])
    if np.isnan(m).all():
        print(f"{roi:5s}  {'—':4s}")
        continue
    pk = int(np.argmax(m))
    peak_lag_TR[roi] = LAGS[pk]
    print(f"{roi:5s}  {n:4d}  {LAGS[pk]:+11d}  {lags_s[pk]:+10.1f}  "
          f"{m[pk]:8.4f}  {m[lag0_idx]:8.4f}")"""
))

# ── Cell 7: figure ────────────────────────────────────────────────────────────
cells.append(code(
"""VISUAL_ROIS = ['LGN', 'V1', 'V2', 'V3V', 'V4', 'V5']
MOTOR_ROIS  = ['IPS', 'FEF', 'SEF', 'PCC']

fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=False)
titles = ['Visual ROIs  (LGN, V1–V5)',
          'Motor / network ROIs  (IPS, FEF, SEF, PCC)']

peak_lag_s = {}
for ax, roi_group, title in zip(axes, [VISUAL_ROIS, MOTOR_ROIS], titles):
    ax.axhline(0, color='k', lw=0.6)
    ax.axvline(0, color='k', lw=0.8, ls='--', alpha=0.35)
    for roi in roi_group:
        m = means[roi]
        s = sems[roi]
        if np.isnan(m).all():
            continue
        c    = ROI_COLORS[roi]
        pk   = int(np.argmax(m))
        pk_s = lags_s[pk]
        peak_lag_s[roi] = pk_s
        ax.plot(lags_s, m, color=c, lw=2.2, label=roi, zorder=3)
        ax.fill_between(lags_s, m - s, m + s, color=c, alpha=0.15, zorder=2)
        ax.scatter([pk_s], [m[pk]], color=c, s=55, zorder=5,
                   edgecolors='white', lw=0.8)
        ax.axvline(pk_s, color=c, lw=0.8, ls=':', alpha=0.6)
    ax.set_xlabel('Lag (s)   [positive = eye precedes BOLD]', fontsize=11)
    ax.set_xticks(np.arange(-20, 22, 4))
    ax.set_title(title, fontsize=11)
    ax.legend(fontsize=9, framealpha=0.85)
    ax.spines[['top', 'right']].set_visible(False)

axes[0].set_ylabel('Mean Pearson r  (±SEM)', fontsize=11)

# FEF-V1 offset: annotate as text in the right panel
if 'FEF' in peak_lag_s and 'V1' in peak_lag_s:
    delta = abs(peak_lag_s['V1'] - peak_lag_s['FEF'])
    axes[1].text(0.03, 0.97,
                 f'FEF peaks {delta:.0f} s earlier than V1',
                 transform=axes[1].transAxes,
                 fontsize=9, va='top', color='#333',
                 bbox=dict(facecolor='white', edgecolor='#ccc', pad=3))

fig.suptitle(f'Eye–BOLD cross-correlation by ROI — Controls (n={len(xcorr_all)})\\n'
             'Nudelman no-BP · 6 mm FWHM · ±SEM across subjects',
             fontsize=11, y=1.01)
fig.tight_layout()

OUT = os.path.join(FIG_DIR, '1_Eye_Brain_lag', 'crosscorr_controls_2panel')
fig.savefig(OUT + '.pdf', bbox_inches='tight')
fig.savefig(OUT + '.png', dpi=150, bbox_inches='tight')
print(f'Saved: {OUT}.png')
plt.show()"""
))

# ── Cell 8: results note ──────────────────────────────────────────────────────
cells.append(md(
"""## Interpretation notes

- Positive peak lag → eye signal precedes BOLD (consistent with eye → cortex direction, but TR=2 s cannot resolve sub-TR neural ordering).
- Peak lag difference FEF vs V1 (annotated in figure): consistent with a motor-to-sensory propagation sequence within the same eye-movement cycle.
- IPS / SEF going negative at later lags = network suppression or post-movement inhibition, not a separate anti-phasic signal.
- These patterns describe the *temporal structure* of eye–BOLD coupling; causal inference requires sub-second imaging or perturbation."""
))

nb = {
    "nbformat": 4,
    "nbformat_minor": 5,
    "metadata": {
        "kernelspec": {"display_name": "Python 3",
                       "language": "python", "name": "python3"},
        "language_info": {"name": "python", "version": "3.10.0"}
    },
    "cells": cells,
}

out = '/home/cemal/Desktop/Opus/OMvar_paper/scratch_1_Eye_Brain_lag.ipynb'
with open(out, 'w') as f:
    json.dump(nb, f, indent=1)
print(f"Written: {out}")
