"""
VPH 2026 — comorbidity heatmap (Step 3)
Y axis (rows):  5 oculomotor metrics + NIHSS total = 6
X axis (cols):  those 6 + ARAT, BOSTON, HVLT, BVMT, POSNER RT = 11
Diagonal cells (same variable) are left blank.
FDR-BH across all off-diagonal tests.
Method: Pearson r, acute stage only.
"""
import pickle, warnings
import numpy as np
import pandas as pd
from scipy import signal as scipy_signal, stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

OUT = '/home/cemal/Desktop/Opus/vph'

# ── 1. LOAD DATA ──────────────────────────────────────────────
with open('/home/cemal/Desktop/Opus/stroke_csv_in_progress.pkl', 'rb') as f:
    df = pickle.load(f)
patients = df[(df.subj_type == 0) & df.inclusion_notes.isna()].copy()
acute = patients[patients.Session == 'acute'].copy()

# ── 2. EXTRACT FEATURES — FIRST RUN ──────────────────────────
def spectral_entropy(x):
    x = x[np.isfinite(x)]
    if len(x) < 8: return np.nan
    _, pxx = scipy_signal.welch(x, nperseg=min(64, len(x)))
    pxx = pxx / (pxx.sum() + 1e-12)
    return float(-np.sum(pxx * np.log(pxx + 1e-12)))

def features_from_run(arr, run_idx=0):
    if not isinstance(arr, np.ndarray) or arr.ndim != 2: return {}
    if arr.shape[1] <= run_idx: return {}
    x = arr[:, run_idx].astype(float)
    v = x[np.isfinite(x)]
    if len(v) < 10: return {}
    ac1 = float(np.corrcoef(v[:-1], v[1:])[0, 1]) if len(v) > 2 else np.nan
    return dict(mean=float(np.mean(v)), var=float(np.var(v, ddof=1)),
                spec_ent=spectral_entropy(v), ac1=ac1)

feat_rows = []
for _, row in acute.iterrows():
    fl = features_from_run(row.get('r_trtrcorr_eye_l'))
    fr = features_from_run(row.get('r_trtrcorr_eye_r'))
    def bil(a, b):
        v = [x for x in [a, b] if np.isfinite(x)]
        return float(np.mean(v)) if v else np.nan
    lr_diff = (fl.get('mean', np.nan) - fr.get('mean', np.nan))
    if not (np.isfinite(fl.get('mean', np.nan)) and np.isfinite(fr.get('mean', np.nan))):
        lr_diff = np.nan
    feat_rows.append({'ID': row['ID'],
        'mean':     bil(fl.get('mean',np.nan),     fr.get('mean',np.nan)),
        'var':      bil(fl.get('var',np.nan),      fr.get('var',np.nan)),
        'spec_ent': bil(fl.get('spec_ent',np.nan), fr.get('spec_ent',np.nan)),
        'ac1':      bil(fl.get('ac1',np.nan),      fr.get('ac1',np.nan)),
        'lr_diff':  lr_diff})

feat_df = pd.DataFrame(feat_rows)
acute = acute.merge(feat_df, on='ID', how='left')

# ── 3. DEFINE ROWS AND COLUMNS ────────────────────────────────
# Y axis: 6 rows
VISUAL_ROWS = [
    ('mean',      'Mean r(t)'),
    ('var',       'Variance'),
    ('spec_ent',  'Spectral entropy'),
    ('ac1',       'AC1'),
    ('lr_diff',   'LR asymmetry'),
    ('nih_total', 'NIHSS total'),
]

# X axis: 11 columns (same 6 + 5 domain outcomes)
DOMAIN_COLS = [
    ('mean',              'Mean r(t)'),
    ('var',               'Variance'),
    ('spec_ent',          'Spectral entropy'),
    ('ac1',               'AC1'),
    ('lr_diff',           'LR asymmetry'),
    ('nih_total',         'NIHSS total'),
    ('raratotal',         'Motor\n(ARAT total)'),
    ('boston_raw',        'Language\n(Boston Naming)'),
    ('hvlt_im',           'Verbal memory\n(HVLT imm.)'),
    ('bvmt_im',           'Visual memory\n(BVMT imm.)'),
    ('pos_rt_disengage',  'Attention\n(Posner diseng. RT)'),
]

N_SHARED = 6   # first 6 cols are the same vars as the rows

# ── 4. PEARSON CORRELATIONS ───────────────────────────────────
rows = []
for vis_col, vis_lbl in VISUAL_ROWS:
    for dom_col, dom_lbl in DOMAIN_COLS:
        # leave diagonal blank
        if vis_col == dom_col:
            rows.append({'vis': vis_col, 'dom': dom_col, 'r': np.nan, 'p': np.nan, 'n': 0, 'diagonal': True})
            continue
        if vis_col not in acute.columns or dom_col not in acute.columns:
            rows.append({'vis': vis_col, 'dom': dom_col, 'r': np.nan, 'p': np.nan, 'n': 0, 'diagonal': False})
            continue
        d = acute[[vis_col, dom_col]].apply(pd.to_numeric, errors='coerce').dropna()
        if len(d) < 10:
            rows.append({'vis': vis_col, 'dom': dom_col, 'r': np.nan, 'p': np.nan, 'n': len(d), 'diagonal': False})
            continue
        r, p = stats.pearsonr(d[vis_col].values, d[dom_col].values)
        rows.append({'vis': vis_col, 'dom': dom_col, 'r': r, 'p': p, 'n': len(d), 'diagonal': False})

corr_df = pd.DataFrame(rows)

# FDR-BH across all off-diagonal tests with valid p
valid = corr_df['p'].notna() & ~corr_df['diagonal']
_, p_fdr, _, _ = multipletests(corr_df.loc[valid, 'p'], method='fdr_bh')
corr_df.loc[valid, 'p_fdr'] = p_fdr

def sig_stars(p):
    if pd.isna(p): return ''
    return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else '†' if p < 0.10 else ''

corr_df['stars'] = corr_df['p_fdr'].apply(sig_stars)

# ── 5. PRINT TABLE ────────────────────────────────────────────
n_tests = valid.sum()
print("\n" + "="*75)
print(f"COMORBIDITY — Pearson r | Acute stage | FDR-BH ({n_tests} off-diagonal tests)")
print("="*75)
print(f"{'Visual measure':<22} {'Domain':<25} {'r':>6} {'p(raw)':>8} {'p(FDR)':>8} {'sig':>4} {'N':>5}")
print("-"*75)
for _, row in corr_df[~corr_df['diagonal']].sort_values('p_fdr').iterrows():
    if pd.isna(row['r']): continue
    flag = ' ◄' if row.get('p_fdr', 1) < 0.05 else ''
    print(f"{row['vis']:<22} {row['dom']:<25} {row['r']:>6.3f} "
          f"{row['p']:>8.4f} {row.get('p_fdr', np.nan):>8.4f} "
          f"{row['stars']:>4}{flag}  n={int(row['n'])}")

corr_df.to_csv(f'{OUT}/vph_comorbidity.csv', index=False, float_format='%.4f')
print(f"\nSaved vph_comorbidity.csv")

# ── 6. HEATMAP ────────────────────────────────────────────────
vis_order  = [c for c, _ in VISUAL_ROWS]
dom_order  = [c for c, _ in DOMAIN_COLS]
vis_labels = [l for _, l in VISUAL_ROWS]
dom_labels = [l for _, l in DOMAIN_COLS]

n_rows = len(vis_order)
n_cols = len(dom_order)

r_mat = np.full((n_rows, n_cols), np.nan)
p_mat = np.full((n_rows, n_cols), np.nan)
n_mat = np.full((n_rows, n_cols), np.nan)
diag_mask = np.zeros((n_rows, n_cols), dtype=bool)

for i, vis in enumerate(vis_order):
    for j, dom in enumerate(dom_order):
        sel = corr_df[(corr_df.vis == vis) & (corr_df.dom == dom)]
        if len(sel):
            r_mat[i, j] = sel.iloc[0]['r']
            p_mat[i, j] = sel.iloc[0].get('p_fdr', np.nan)
            n_mat[i, j] = sel.iloc[0]['n']
            if sel.iloc[0]['diagonal']:
                diag_mask[i, j] = True

fig, ax = plt.subplots(figsize=(14, 5.5), constrained_layout=True)

# Plot non-diagonal cells; show diagonal as grey
r_plot = r_mat.copy()
im = ax.imshow(r_plot, cmap='RdBu_r', vmin=-0.6, vmax=0.6, aspect='auto')
plt.colorbar(im, ax=ax, label='Pearson r', fraction=0.025)

# Grey out diagonal cells
for i in range(n_rows):
    for j in range(n_cols):
        if diag_mask[i, j]:
            ax.add_patch(plt.Rectangle((j-0.5, i-0.5), 1, 1,
                                       color='#cccccc', zorder=2))
            ax.text(j, i, '—', ha='center', va='center', fontsize=10,
                    color='#888888', zorder=3)

ax.set_xticks(range(n_cols))
ax.set_yticks(range(n_rows))
ax.set_xticklabels(dom_labels, fontsize=9)
ax.set_yticklabels(vis_labels, fontsize=10)

# Cell annotations: r value + FDR stars + N (off-diagonal only)
for i in range(n_rows):
    for j in range(n_cols):
        if diag_mask[i, j]: continue
        r_val = r_mat[i, j]
        if np.isnan(r_val): continue
        star  = sig_stars(p_mat[i, j])
        n     = int(n_mat[i, j])
        color = 'white' if abs(r_val) > 0.35 else 'black'
        weight = 'bold' if star else 'normal'
        ax.text(j, i, f'{r_val:.2f}{star}\nn={n}',
                ha='center', va='center', fontsize=8,
                color=color, fontweight=weight, zorder=4)

# Vertical separator between shared-variable block and domain outcome block
ax.axvline(N_SHARED - 0.5, color='white', linewidth=2.5, zorder=5)
ax.axvline(N_SHARED - 0.5, color='#444444', linewidth=1.5, linestyle='--', zorder=6)

# Column group labels on top
ax.text(2.5, -1.1, 'Oculomotor metrics / NIHSS', transform=ax.transData,
        fontsize=8.5, color='#333333', va='center', ha='center', style='italic')
ax.text(8.0, -1.1, 'Domain outcomes', transform=ax.transData,
        fontsize=8.5, color='#333333', va='center', ha='center', style='italic')

ax.set_title('Comorbidity — acute oculomotor metrics × functional domain outcomes\n'
             f'Pearson r  |  FDR-BH corrected ({n_tests} off-diagonal tests)  |  Acute stage only',
             fontsize=11, fontweight='bold')
ax.spines[['top','right','bottom','left']].set_visible(False)

plt.savefig(f'{OUT}/vph_comorbidity_heatmap.png', dpi=200, bbox_inches='tight', facecolor='white')
plt.close()
print("Saved vph_comorbidity_heatmap.png")
print("\nDone.")
