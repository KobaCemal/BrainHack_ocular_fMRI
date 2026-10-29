#!/usr/bin/env python3
"""
Surface-space merged GLM: Baltimore (SC) + Canada (SC).
Both datasets have micapipe fsLR-32k cleaned surface timeseries.
No alignment needed — same 64984-vertex space.
GLM: eye_hrf ~ surface timeseries (OLS, intercept only).
Group model: intercept + site (Canada=1, Baltimore=0).
"""
import os, pickle, re
import numpy as np
import nibabel as nib
from scipy.signal import convolve
from scipy.stats import gamma
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from nilearn import plotting, surface
from tqdm import tqdm

BALT_DIR      = '/mnt/compneuro/MULTIBOOT/blindness_gradients/datasets/baltimore'
CANA_DIR      = '/mnt/compneuro/MULTIBOOT/blindness_gradients/datasets/canada'
BALT_MICAPIPE = os.path.join(BALT_DIR, 'derivatives/micapipe_v0.2.0')
CANA_MICAPIPE = os.path.join(CANA_DIR, 'derivatives/micapipe_v0.2.0')
BALT_TRC_PKL  = os.path.join(BALT_DIR, 'derivatives/deepmreye/trtrcorr_baltimore_sc.pkl')
CANA_TRC_PKL  = os.path.join(CANA_DIR, 'derivatives/deepmreye/trtrcorr_canada_sc.pkl')
GLASSER_R     = '/mnt/compneuro/MULTIBOOT/blindness_gradients/source/parcellations/glasser_in_fslr32_r.gii'
GLASSER_L     = '/mnt/compneuro/MULTIBOOT/blindness_gradients/source/parcellations/glasser_in_fslr32_l.gii'
OUT_DIR       = '/tmp/merged_surface'
os.makedirs(OUT_DIR, exist_ok=True)

TR_BALT = 2.0
TR_CANA = 2.2
OUTLIERS_BALT = {'sub-029', 'sub-051'}


def spm_hrf(tr, duration=32.0):
    t = np.arange(0, duration, tr)
    h = gamma.pdf(t, 6, scale=1.0) - gamma.pdf(t, 16, scale=1.0) / 6
    return h / np.abs(h).max()


def hrf_convolve(ts, hrf, n_trs):
    c   = np.nan_to_num(ts, nan=0.0)
    out = convolve(c, hrf, mode='full')[:n_trs]
    return out - out.mean()


def run_surface_glm(surf_ts, eye_ts, tr):
    n_trs  = surf_ts.shape[0]
    hrf    = spm_hrf(tr)
    eye_r  = hrf_convolve(eye_ts, hrf, n_trs)
    # z-score so beta is comparable across subjects regardless of eye-signal scale
    sd = eye_r.std()
    if sd < 1e-10:
        return np.zeros(surf_ts.shape[1])
    eye_r  = (eye_r - eye_r.mean()) / sd
    X      = np.column_stack([eye_r, np.ones(n_trs)])
    betas, _, _, _ = np.linalg.lstsq(X, surf_ts, rcond=None)
    return betas[0]   # (64984,) cope


# Glasser V1 mask
r_labels = nib.load(GLASSER_R).darrays[0].data  # LH parcels 1-180
l_labels = nib.load(GLASSER_L).darrays[0].data  # RH parcels 181-360
glasser  = np.concatenate([r_labels, l_labels])  # (64984,)
v1_mask  = np.isin(glasser, [1, 181])
print(f'V1 vertices: {v1_mask.sum()}')

# Load PKLs
with open(BALT_TRC_PKL, 'rb') as f:
    balt_trc = pickle.load(f)
with open(CANA_TRC_PKL, 'rb') as f:
    cana_trc = pickle.load(f)

records = []

# Baltimore — strict session matching
seen_balt = set()
for run_id in sorted(balt_trc.keys()):
    d   = balt_trc[run_id]
    sub = d['sub']
    if sub in OUTLIERS_BALT or sub in seen_balt:
        continue
    m   = re.search(r'(ses-\d+)', run_id)
    ses = m.group(1) if m else None
    if not ses:
        continue
    surf = os.path.join(BALT_MICAPIPE, sub, ses, 'func',
                        'desc-se_task-rest_bold', 'surf',
                        f'{sub}_{ses}_surf-fsLR-32k_desc-timeseries_clean.shape.gii')
    if not os.path.exists(surf):
        continue
    seen_balt.add(sub)
    records.append({'sub': sub, 'site': 0, 'tr': TR_BALT,
                    'surf_path': surf, 'eye_ts': d['r_timeseries']})

n_balt = len(records)
print(f'Baltimore: {n_balt}')

# Canada — no session in path
seen_cana = set()
for run_id in sorted(cana_trc.keys()):
    d   = cana_trc[run_id]
    sub = d['sub']
    if sub in seen_cana:
        continue
    surf = os.path.join(CANA_MICAPIPE, sub, 'func',
                        'desc-se_task-rest_bold', 'surf',
                        f'{sub}_surf-fsLR-32k_desc-timeseries_clean.shape.gii')
    if not os.path.exists(surf):
        continue
    seen_cana.add(sub)
    records.append({'sub': sub, 'site': 1, 'tr': TR_CANA,
                    'surf_path': surf, 'eye_ts': d['r_timeseries']})

print(f'Canada: {len(records) - n_balt}')
print(f'Total: {len(records)}')

# Run GLMs
cope_maps = []
v1_betas  = []
sites     = []
subs      = []

for rec in tqdm(records):
    try:
        surf_ts_full = nib.load(rec['surf_path']).darrays[0].data  # (n_trs, 64984)
        surf_ts = surf_ts_full[1:-1, :]                            # trim first/last TR
        n_trs   = surf_ts.shape[0]
        eye_ts  = rec['eye_ts'][:n_trs]
        cope    = run_surface_glm(surf_ts, eye_ts, rec['tr'])
        cope_maps.append(cope)
        v1_betas.append(float(cope[v1_mask].mean()))
        sites.append(rec['site'])
        subs.append(rec['sub'])
        print(f"  {rec['sub']} (site={rec['site']}): V1 beta={v1_betas[-1]:.4f}")
    except Exception as e:
        print(f"  ERROR {rec['sub']}: {e}")

cope_maps = np.array(cope_maps)
sites     = np.array(sites)
v1_betas  = np.array(v1_betas)
print(f'\nCope maps: {cope_maps.shape}')
print(f'V1 beta Baltimore: {v1_betas[sites==0].mean():.4f} ± {v1_betas[sites==0].std():.4f}')
print(f'V1 beta Canada:    {v1_betas[sites==1].mean():.4f} ± {v1_betas[sites==1].std():.4f}')

np.save(os.path.join(OUT_DIR, 'cope_maps.npy'), cope_maps)
np.save(os.path.join(OUT_DIR, 'sites.npy'), sites)
np.save(os.path.join(OUT_DIR, 'v1_betas.npy'), v1_betas)

# Group OLS: intercept + site
n_sub    = cope_maps.shape[0]
X_grp    = np.column_stack([np.ones(n_sub), sites])
betas_g, _, _, _ = np.linalg.lstsq(X_grp, cope_maps, rcond=None)
resid    = cope_maps - X_grp @ betas_g
sigma2   = (resid ** 2).sum(axis=0) / (n_sub - 2)
XtXinv   = np.linalg.inv(X_grp.T @ X_grp)
t_int    = betas_g[0] / np.sqrt(XtXinv[0, 0] * sigma2 + 1e-12)
t_site   = betas_g[1] / np.sqrt(XtXinv[1, 1] * sigma2 + 1e-12)

np.save(os.path.join(OUT_DIR, 't_intercept.npy'), t_int)
np.save(os.path.join(OUT_DIR, 't_site.npy'), t_site)
print(f'V1 t_intercept: mean={t_int[v1_mask].mean():.3f}, max={t_int[v1_mask].max():.3f}')

# Save surface t-maps as GIFTIs (LH/RH split)
def save_gifti(data, path):
    arr  = nib.gifti.GiftiDataArray(data.astype(np.float32))
    img  = nib.gifti.GiftiImage(darrays=[arr])
    nib.save(img, path)

save_gifti(t_int[:32492],  os.path.join(OUT_DIR, 't_intercept_LH.func.gii'))
save_gifti(t_int[32492:],  os.path.join(OUT_DIR, 't_intercept_RH.func.gii'))
save_gifti(t_site[:32492], os.path.join(OUT_DIR, 't_site_LH.func.gii'))
save_gifti(t_site[32492:], os.path.join(OUT_DIR, 't_site_RH.func.gii'))
print('Saved GIFTI t-maps')

# Bar plot: V1 beta per site (z-scored betas, now comparable)
fig, ax = plt.subplots(figsize=(5, 4))
np.random.seed(42)
for i, (label, mask, col) in enumerate([
        ('Baltimore', sites==0, 'steelblue'),
        ('Canada',    sites==1, 'coral')]):
    vals = v1_betas[mask]
    ax.bar(i, vals.mean(), yerr=vals.std()/np.sqrt(len(vals)), capsize=5,
           color=col, alpha=0.8)
    ax.scatter(np.full(len(vals), i) + np.random.uniform(-0.15, 0.15, len(vals)),
               vals, c='k', s=15, zorder=3)
ax.set_xticks([0, 1])
ax.set_xticklabels([f'Baltimore\n(n={(sites==0).sum()})',
                    f'Canada\n(n={(sites==1).sum()})'])
ax.set_ylabel('V1 beta (z-scored eye regressor)')
ax.set_title('Eye–V1 coupling (surface GLM)')
ax.axhline(0, c='k', lw=0.7, ls='--')
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, 'v1_beta_bar.png'), dpi=150)
plt.close()
print('Saved v1_beta_bar.png')

# Histogram of V1 betas
fig, ax = plt.subplots(figsize=(6, 4))
ax.hist(v1_betas[sites==0], bins=12, alpha=0.6, color='steelblue', label='Baltimore')
ax.hist(v1_betas[sites==1], bins=6,  alpha=0.6, color='coral',     label='Canada')
ax.axvline(0, c='k', lw=0.8, ls='--')
ax.set_xlabel('V1 beta (z-scored eye regressor)')
ax.set_ylabel('Count')
ax.legend()
ax.set_title('V1 beta distribution per site')
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, 'v1_beta_hist.png'), dpi=150)
plt.close()
print('Saved v1_beta_hist.png')
print('Done.')
