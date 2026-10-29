#!/usr/bin/env python3
"""
Compute r_trtrcorr for Canada sighted controls.
Confounds: 6 motion params from micapipe desc-se.1D + derivatives + FD
           + CSF, WM, GM tissue signals (16 cols total).
"""
import os, glob, pickle, re
import numpy as np
import pandas as pd
from nilearn import signal as nl_signal
from tqdm import tqdm

BIDS_DIR   = '/mnt/compneuro/MULTIBOOT/blindness_gradients/datasets/canada'
DME_DIR    = os.path.join(BIDS_DIR, 'derivatives/deepmreye')
MICAPIPE   = os.path.join(BIDS_DIR, 'derivatives/micapipe_v0.2.0')
OUTPUT_PKL = os.path.join(DME_DIR, 'trtrcorr_canada_sc.pkl')

EYE_R_SLICE = slice(0, 23)
EYE_L_SLICE = slice(24, 47)


def build_confounds(sub, n_trs):
    func_dir = os.path.join(MICAPIPE, sub, 'func',
                            'desc-se_task-rest_bold', 'volumetric')
    mot_path = os.path.join(func_dir, f'{sub}_space-func_desc-se.1D')
    params   = np.loadtxt(mot_path)[:n_trs]   # (n_trs, 6)
    derivs   = np.vstack([np.zeros(6), np.diff(params, axis=0)])
    rot_d, tx_d = derivs[:, :3], derivs[:, 3:]
    fd       = np.abs(rot_d * 50).sum(1) + np.abs(tx_d).sum(1)
    fd[0]    = 0.0
    conf     = np.hstack([params, derivs, fd.reshape(-1, 1)])   # 13 cols

    for tissue in ('CSF', 'WM', 'GM'):
        txt = os.path.join(func_dir, f'{sub}_space-func_desc-se_pve_{tissue}.txt')
        if os.path.exists(txt):
            ts = np.loadtxt(txt)[:n_trs].reshape(-1, 1)
            if len(ts) < n_trs:
                ts = np.pad(ts, ((0, n_trs - len(ts)), (0, 0)))
            conf = np.hstack([conf, ts])

    return conf   # 16 cols


def compute_trtrcorr(data4d):
    n_trs    = data4d.shape[-1]
    r_series = np.full(n_trs, np.nan)
    flat     = data4d.reshape(-1, n_trs)
    for t in range(1, n_trs):
        v1, v2 = flat[:, t - 1], flat[:, t]
        valid  = np.isfinite(v1) & np.isfinite(v2)
        if valid.sum() > 10:
            r_series[t] = np.corrcoef(v1[valid], v2[valid])[0, 1]
    return r_series


def residualize_and_clamp(om_slice, conf):
    signal      = om_slice.copy()
    nm          = ~np.isfinite(signal); signal[nm] = 0.0
    conf        = np.nan_to_num(conf.astype(np.float64), nan=0.0)
    cleaned     = nl_signal.clean(signal.reshape(-1, 1), confounds=conf,
                                  detrend=False, standardize=False).flatten()
    cleaned[nm] = 0.0
    return np.maximum(cleaned, 0.0)


def two_point_ma(arr):
    return np.array([0.5 * (arr[i] + arr[i + 1]) for i in range(len(arr) - 1)])


p_files = (sorted(glob.glob(os.path.join(DME_DIR, 'sub-*', 'mask_*_preproc.p'))) +
           sorted(glob.glob(os.path.join(DME_DIR, 'sub-*', 'mask_*_mni.p'))))
print(f'Found {len(p_files)} .p files')

results = {}
if os.path.exists(OUTPUT_PKL):
    with open(OUTPUT_PKL, 'rb') as f:
        results = pickle.load(f)
    print(f'Resuming — {len(results)} done')

n_done = n_skip = n_err = 0

for p_path in tqdm(p_files):
    fname  = os.path.basename(p_path)
    run_id = re.sub(r'\.p$', '', fname.replace('mask_', ''))

    if run_id in results:
        n_skip += 1; continue

    m = re.search(r'(sub-s(?:ccb|clb)\d+)', fname)
    if not m:
        print(f'Cannot parse sub from {fname}'); n_err += 1; continue
    sub = m.group(1)

    try:
        with open(p_path, 'rb') as f:
            data4d = np.array(pickle.load(f), dtype=np.float32)

        n_trs = data4d.shape[3]
        conf  = build_confounds(sub, n_trs)
        ct    = conf[1:]

        r_both  = compute_trtrcorr(data4d)
        r_eye_r = compute_trtrcorr(data4d[EYE_R_SLICE])
        r_eye_l = compute_trtrcorr(data4d[EYE_L_SLICE])

        om_both_rc  = residualize_and_clamp(1.0 - r_both[1:],  ct)
        om_eye_r_rc = residualize_and_clamp(1.0 - r_eye_r[1:], ct)
        om_eye_l_rc = residualize_and_clamp(1.0 - r_eye_l[1:], ct)

        results[run_id] = {
            'sub':          sub,
            'n_trs':        n_trs,
            'r_timeseries': two_point_ma(om_both_rc),
            'mean_r':       float(np.nanmean(two_point_ma(om_both_rc))),
            'var_r':        float(np.nanvar(two_point_ma(om_both_rc))),
            'r_eye_r':      two_point_ma(om_eye_r_rc),
            'r_eye_l':      two_point_ma(om_eye_l_rc),
        }
        n_done += 1
        if n_done % 5 == 0:
            with open(OUTPUT_PKL, 'wb') as f:
                pickle.dump(results, f)

    except Exception as e:
        print(f'ERROR {run_id}: {e}'); n_err += 1

print(f'\nDone: {n_done} | Skipped: {n_skip} | Errors: {n_err}')
with open(OUTPUT_PKL, 'wb') as f:
    pickle.dump(results, f)
print(f'Saved → {OUTPUT_PKL}')
