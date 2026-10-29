#!/usr/bin/env python3
"""DeepMReye eye mask extraction for Canada sighted controls (18 unique subjects).
Uses micapipe preproc BOLD directly — no motion correction step needed.
Symlinks preproc NIfTIs into deepmreye dir so .p files land there.
"""
import os, glob
from multiprocessing import Pool
from tqdm import tqdm

BIDS_DIR  = '/mnt/compneuro/MULTIBOOT/blindness_gradients/datasets/canada'
MICAPIPE  = os.path.join(BIDS_DIR, 'derivatives/micapipe_v0.2.0')
DME_DIR   = os.path.join(BIDS_DIR, 'derivatives/deepmreye')
N_WORKERS = 6

SIGHTED = [
    'sub-sccb01', 'sub-sccb02', 'sub-sccb03', 'sub-sccb04', 'sub-sccb05',
    'sub-sccb06', 'sub-sccb07', 'sub-sccb08', 'sub-sccb09', 'sub-sccb10',
    'sub-sccb11', 'sub-sccb12', 'sub-sccb13', 'sub-sccb14',
    'sub-sclb01', 'sub-sclb04', 'sub-sclb06', 'sub-sclb11',
]

# ── Symlink preproc NIfTIs into deepmreye dir ─────────────────────────────────
all_links = []
for sub in SIGHTED:
    src = os.path.join(MICAPIPE, sub, 'func', 'desc-se_task-rest_bold',
                       'volumetric', f'{sub}_space-func_desc-se_preproc.nii.gz')
    if not os.path.exists(src):
        print(f'Missing preproc: {sub}')
        continue
    out_dir = os.path.join(DME_DIR, sub)
    os.makedirs(out_dir, exist_ok=True)
    link = os.path.join(out_dir, f'{sub}_task-rest_bold_preproc.nii.gz')
    if not os.path.exists(link):
        os.symlink(src, link)
    all_links.append(link)

print(f'Found {len(all_links)} subjects\n')


def init_worker():
    global eyemask_small, eyemask_big, dme_template, mask, x_edges, y_edges, z_edges
    from deepmreye import preprocess
    (eyemask_small, eyemask_big, dme_template,
     mask, x_edges, y_edges, z_edges) = preprocess.get_masks()


def process_run(fp_func):
    from deepmreye import preprocess
    base_name   = os.path.basename(fp_func).replace('.nii.gz', '')
    output_mask = os.path.join(os.path.dirname(fp_func), f'mask_{base_name}.p')
    if os.path.exists(output_mask):
        return fp_func, 'skip'
    try:
        preprocess.run_participant(
            fp_func, dme_template, eyemask_big, eyemask_small,
            x_edges, y_edges, z_edges,
            transforms=['Affine', 'Affine', 'SyNAggro'])
        return fp_func, 'done'
    except Exception as e:
        return fp_func, f'error: {e}'


n_done = n_skip = n_err = 0
with Pool(N_WORKERS, initializer=init_worker) as pool:
    for fp, status in tqdm(pool.imap_unordered(process_run, all_links), total=len(all_links)):
        if status == 'done':
            print(f'  Done: {os.path.basename(fp)}'); n_done += 1
        elif status == 'skip':
            n_skip += 1
        else:
            print(f'  {status} — {os.path.basename(fp)}'); n_err += 1

print(f'\nDone: {n_done} | Skipped: {n_skip} | Errors: {n_err}')
