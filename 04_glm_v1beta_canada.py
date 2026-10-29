#!/usr/bin/env python3
"""
Nudelman no-BP GLM for Canada sighted controls.
Eye predictor: r_trtrcorr HRF-convolved.
BOLD input: micapipe preproc NIfTI (already motion-corrected).
"""
import os, pickle, glob, re, tempfile
import numpy as np
import pandas as pd
import nibabel as nib
import xml.etree.ElementTree as ET
from scipy.signal import convolve
from scipy.stats import gamma
from scipy.ndimage import affine_transform
from nilearn.glm.first_level import FirstLevelModel, make_first_level_design_matrix
from nilearn.glm.second_level import SecondLevelModel
from tqdm import tqdm

BIDS_DIR     = '/mnt/compneuro/MULTIBOOT/blindness_gradients/datasets/canada'
DME_DIR      = os.path.join(BIDS_DIR, 'derivatives/deepmreye')
MICAPIPE     = os.path.join(BIDS_DIR, 'derivatives/micapipe_v0.2.0')
os.makedirs(DME_DIR, exist_ok=True)
OUT_DIR      = os.path.join(BIDS_DIR, 'derivatives/glm_nudelman')
TRTRCORR_PKL = os.path.join(DME_DIR, 'trtrcorr_canada_sc.pkl')
OUTPUT_PKL   = os.path.join(OUT_DIR, 'v1beta_canada_sc.pkl')
ATLAS_PATH   = '/mnt/compneuro/fsl/data/atlases/Juelich/Juelich-maxprob-thr25-2mm.nii.gz'
FSL_ATLASES  = '/mnt/compneuro/fsl/data/atlases'

TR      = 2.2
FWHM_MM = 6.0
os.makedirs(OUT_DIR, exist_ok=True)

tree   = ET.parse(os.path.join(FSL_ATLASES, 'Juelich.xml'))
labels = [(int(el.get('index')), el.text.strip()) for el in tree.findall('.//label')]
V1_IDX = [i + 1 for i, name in labels if 'V1' in name]


def spm_hrf(tr=TR, duration=32.0):
    t = np.arange(0, duration, tr)
    h = gamma.pdf(t, 6, scale=1.0) - gamma.pdf(t, 16, scale=1.0) / 6
    return h / np.abs(h).max()

hrf = spm_hrf()

def hrf_convolve(ts, n_trs):
    c   = np.nan_to_num(ts, nan=0.0)
    out = convolve(c, hrf, mode='full')[:n_trs]
    return out - out.mean()

NUDELMAN_COLS = [
    'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z',
    'trans_x_derivative1', 'trans_y_derivative1', 'trans_z_derivative1',
    'rot_x_derivative1', 'rot_y_derivative1', 'rot_z_derivative1',
    'framewise_displacement', 'csf', 'white_matter', 'global_signal',
]

def build_confounds_df(sub, n_trs):
    func_dir = os.path.join(MICAPIPE, sub, 'func',
                            'desc-se_task-rest_bold', 'volumetric')
    mot_path = os.path.join(func_dir, f'{sub}_space-func_desc-se.1D')
    params   = np.loadtxt(mot_path)[:n_trs]
    derivs   = np.vstack([np.zeros(6), np.diff(params, axis=0)])
    rot_d, tx_d = derivs[:, :3], derivs[:, 3:]
    fd       = np.abs(rot_d * 50).sum(1) + np.abs(tx_d).sum(1)
    fd[0]    = 0.0

    tissue = {}
    for name, key in [('CSF', 'csf'), ('WM', 'white_matter'), ('GM', 'global_signal')]:
        txt = os.path.join(func_dir, f'{sub}_space-func_desc-se_pve_{name}.txt')
        ts  = np.loadtxt(txt)[:n_trs] if os.path.exists(txt) else np.zeros(n_trs)
        tissue[key] = ts

    df = pd.DataFrame({
        'rot_x': params[:, 0], 'rot_y': params[:, 1], 'rot_z': params[:, 2],
        'trans_x': params[:, 3], 'trans_y': params[:, 4], 'trans_z': params[:, 5],
        'rot_x_derivative1': derivs[:, 0], 'rot_y_derivative1': derivs[:, 1],
        'rot_z_derivative1': derivs[:, 2],
        'trans_x_derivative1': derivs[:, 3], 'trans_y_derivative1': derivs[:, 4],
        'trans_z_derivative1': derivs[:, 5],
        'framewise_displacement': fd,
        'csf': tissue['csf'], 'white_matter': tissue['white_matter'],
        'global_signal': tissue['global_signal'],
    })
    return df[NUDELMAN_COLS]


def warp_v1_to_func(sub, ref_img_path, tmp_dir):
    atlas_img = nib.load(ATLAS_PATH)
    ref_img   = nib.load(ref_img_path)
    M = np.linalg.inv(atlas_img.affine) @ ref_img.affine
    resampled = affine_transform(
        atlas_img.get_fdata(dtype=np.float32), M[:3, :3],
        offset=M[:3, 3], output_shape=ref_img.shape[:3], order=0, cval=0.0
    ).astype(np.int32)
    out_path = os.path.join(tmp_dir, f'{sub}_v1_func.nii.gz')
    nib.save(nib.Nifti1Image(resampled, ref_img.affine), out_path)
    return out_path


with open(TRTRCORR_PKL, 'rb') as f:
    trtrcorr = pickle.load(f)
print(f'Loaded trtrcorr: {len(trtrcorr)} subjects')

results = {}
if os.path.exists(OUTPUT_PKL):
    with open(OUTPUT_PKL, 'rb') as f:
        results = pickle.load(f)
    print(f'Resuming — {len(results)} done')

n_done = n_skip = n_err = 0

with tempfile.TemporaryDirectory() as tmp_dir:
    for run_id, d in tqdm(trtrcorr.items()):
        sub = d['sub']

        if sub in results:
            n_skip += 1; continue

        func_dir   = os.path.join(MICAPIPE, sub, 'func',
                                  'desc-se_task-rest_bold', 'volumetric')
        preproc    = os.path.join(func_dir, f'{sub}_space-func_desc-se_preproc.nii.gz')
        mni_bold   = os.path.join(DME_DIR, sub, f'{sub}_task-rest_bold_mni.nii.gz')
        bold_path  = preproc if os.path.exists(preproc) and '_mni' not in run_id else mni_bold
        cope_out   = os.path.join(OUT_DIR, f'{sub}_cope.nii.gz')

        if not os.path.exists(bold_path):
            print(f'Missing BOLD: {sub}'); n_err += 1; continue

        try:
            img   = nib.load(bold_path)
            data  = img.get_fdata(dtype=np.float32)[:, :, :, 1:-1]
            n_trs = data.shape[3]

            # explicit brain mask: voxels with non-zero variance across time
            mask_arr = (data.std(axis=3) > 0).astype(np.int8)
            mask_img = nib.Nifti1Image(mask_arr, img.affine)

            eye_ts  = d['r_timeseries'][:n_trs]
            conf_df = build_confounds_df(sub, n_trs + 2)
            conf_df = conf_df.iloc[1:-1].reset_index(drop=True).fillna(0.0)
            conf_df['eye_regressor'] = hrf_convolve(eye_ts, n_trs)

            dm = make_first_level_design_matrix(
                np.arange(n_trs) * TR,
                hrf_model=None, drift_model='cosine',
                add_regs=conf_df.values,
                add_reg_names=list(conf_df.columns))

            fmodel = FirstLevelModel(
                t_r=TR, noise_model='ar1', standardize=False,
                smoothing_fwhm=FWHM_MM, mask_img=mask_img,
                minimize_memory=True, verbose=0)
            fmodel.fit(nib.Nifti1Image(data, img.affine, img.header),
                       design_matrices=dm)
            cope_img = fmodel.compute_contrast('eye_regressor', output_type='effect_size')
            cope_img.to_filename(cope_out)

            v1_func  = warp_v1_to_func(sub, bold_path, tmp_dir)
            v1_arr   = nib.load(v1_func).get_fdata(dtype=np.float32)
            v1_bool  = np.isin(v1_arr, V1_IDX)
            cope_data = cope_img.get_fdata(dtype=np.float32)
            v1_beta   = float(cope_data[v1_bool].mean()) if v1_bool.sum() > 0 else np.nan

            results[sub] = {'sub': sub, 'v1_beta': v1_beta, 'cope_path': cope_out}
            n_done += 1
            print(f'  {sub}: V1 beta={v1_beta:.4f}')

            if n_done % 5 == 0:
                with open(OUTPUT_PKL, 'wb') as f:
                    pickle.dump(results, f)

        except Exception as e:
            print(f'  ERROR {sub}: {e}'); n_err += 1

print(f'\nDone: {n_done} | Skipped: {n_skip} | Errors: {n_err}')
with open(OUTPUT_PKL, 'wb') as f:
    pickle.dump(results, f)

cope_imgs = [nib.load(v['cope_path']) for v in results.values() if os.path.exists(v['cope_path'])]
if cope_imgs:
    dm_grp = pd.DataFrame([1] * len(cope_imgs), columns=['intercept'])
    m = SecondLevelModel(n_jobs=4)
    m.fit(cope_imgs, design_matrix=dm_grp)
    m.compute_contrast(output_type='stat').to_filename(os.path.join(OUT_DIR, 'group_tmap_canada_sc.nii.gz'))
    m.compute_contrast(output_type='z_score').to_filename(os.path.join(OUT_DIR, 'group_zmap_canada_sc.nii.gz'))
    print(f'Group t-map saved → {OUT_DIR}/group_tmap_canada_sc.nii.gz')

print(f'\nV1 beta PKL saved → {OUTPUT_PKL}')
