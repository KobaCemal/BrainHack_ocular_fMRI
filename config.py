"""
Shared configuration for OMvar paper notebooks.
Update here → change propagates to all notebooks that import this file.

Scratch notebooks: use scratch_<topic>.ipynb — never treated as final modules.
"""
import os
import numpy as np
import pandas as pd

# ── Paths ─────────────────────────────────────────────────────────────────────
PKL       = '/mnt/compneuro/deepmreye_finetuning/stroke_csv_in_progress.pkl'
FMRIPREP  = ('/mnt/compneuro/MULTIBOOT/net/ascratch/people/plgkoba'
             '/stroke_BIDS_first_sessions/derivatives/fmriprep')
COPE_BASE = '/mnt/compneuro/deepmreye_finetuning/compare4'

NB_DIR  = '/mnt/compneuro/comparison_of_methods/OMvar_paper'
FIG_DIR = os.path.join(NB_DIR, 'figures')
INT_DIR = os.path.join(NB_DIR, 'intermediates')

# Beta maps: Nudelman no-BP pipeline (see backup/3_SeedBasedGLM_lesionsize.ipynb)
# Per-subject cope files: {NUDELMAN_BETA_DIR}/{sub_id}_cope_fe.nii.gz
NUDELMAN_TMAP     = os.path.join(COPE_BASE, 'nudelman_group_tmap.nii.gz')
NUDELMAN_BETA_DIR = os.path.join(COPE_BASE, 'nudelman_noop')

# ── Acquisition ───────────────────────────────────────────────────────────────
TR      = 2.0
FWHM_MM = 6.0

# ── Confound columns: Nudelman no-BP (16 params) ─────────────────────────────
NUDELMAN_COLS = [
    'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z',
    'trans_x_derivative1', 'trans_y_derivative1', 'trans_z_derivative1',
    'rot_x_derivative1', 'rot_y_derivative1', 'rot_z_derivative1',
    'framewise_displacement', 'csf', 'white_matter', 'global_signal',
]

# ── Group definitions ─────────────────────────────────────────────────────────
GROUP_COLORS = {
    'Controls': '#2166ac',
    'Q1':       '#7ec87e',
    'Q2Q3':     '#e8c45e',
    'Q4':       '#d62728',
}

def get_groups(meta):
    """
    Returns (groups_dict, q25, q75).
    groups_dict: {name: DataFrame of rows from meta}

    Controls : ses-control, subj_type == 1
    Q1       : acute, lesion_volume_mm3 <= 25th pct  (smallest)
    Q2Q3     : acute, 25th < lesion_volume_mm3 <= 75th pct
    Q4       : acute, lesion_volume_mm3 > 75th pct   (largest)

    Exactly as computed in 0_ClinicalSummary.ipynb.
    """
    controls = meta[
        (meta['Session'] == 'control') & (meta['subj_type'] == 1.0)
    ].copy()

    acute = meta[meta['Session'] == 'acute'].copy()
    acute['lesion_volume_mm3'] = pd.to_numeric(
        acute['lesion_volume_mm3'], errors='coerce'
    )
    has_vol = acute['lesion_volume_mm3'].notna()
    q25, q75 = acute.loc[has_vol, 'lesion_volume_mm3'].quantile([0.25, 0.75])

    q1_ids   = acute.loc[has_vol & (acute['lesion_volume_mm3'] <= q25),
                         'ID'].unique()
    q2q3_ids = acute.loc[has_vol &
                         (acute['lesion_volume_mm3'] > q25) &
                         (acute['lesion_volume_mm3'] <= q75),
                         'ID'].unique()
    q4_ids   = acute.loc[has_vol & (acute['lesion_volume_mm3'] > q75),
                         'ID'].unique()

    def acute_rows(ids):
        return meta[
            (meta['Session'] == 'acute') & meta['ID'].isin(ids)
        ].copy()

    groups = {
        'Controls': controls,
        'Q1':       acute_rows(q1_ids),
        'Q2Q3':     acute_rows(q2q3_ids),
        'Q4':       acute_rows(q4_ids),
    }
    return groups, float(q25), float(q75)

# ── ROI definitions ───────────────────────────────────────────────────────────
# Juelich maxprob-thr0-1mm atlas labels
JUELICH_LABELS = {
    'LGN':  ['GM Lateral geniculate body'],
    'V1':   ['GM Visual cortex V1 BA17'],
    'V2':   ['GM Visual cortex V2 BA18'],
    'V3V':  ['GM Visual cortex V3V'],
    'V4':   ['GM Visual cortex V4'],
    'V5':   ['GM Visual cortex V5'],
    'IPS':  ['GM Anterior intra-parietal sulcus hIP1',
             'GM Anterior intra-parietal sulcus hIP2',
             'GM Anterior intra-parietal sulcus hIP3'],
}

# Sphere ROIs: FEF bilateral ±28,−6,+50 mm r=8mm (Amiez & Petrides 2009)
#              SEF midline 0,+4,+52 mm r=8mm
#              PCC midline 0,−52,+26 mm r=8mm
SPHERE_ROIS = {
    'FEF': {'centers': [(28, -6, 50), (-28, -6, 50)], 'radius': 8},
    'SEF': {'centers': [(0,  4, 52)],                  'radius': 8},
    'PCC': {'centers': [(0, -52, 26)],                 'radius': 8},
}

ROI_ORDER = ['LGN', 'V1', 'V2', 'V3V', 'V4', 'V5', 'IPS', 'FEF', 'SEF', 'PCC']

ROI_COLORS = {
    'LGN':  '#6a0dad',
    'V1':   '#003f88',
    'V2':   '#0070c0',
    'V3V':  '#00b0f0',
    'V4':   '#00b050',
    'V5':   '#70ad47',
    'IPS':  '#c55a11',
    'FEF':  '#d62728',
    'SEF':  '#ff7f0e',
    'PCC':  '#7b2d8b',
}

# ── Shared helpers ────────────────────────────────────────────────────────────
def parse_runs(row):
    nums = [int(x) for x in str(row.get('Run', '')).split('-')
            if x.strip().isdigit()]
    return list(enumerate(nums))
