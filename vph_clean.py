"""
VPH 2026 — clean hybrid analysis
Features extracted from FIRST RUN ONLY at each session:
  mean, variance, spectral entropy, AC1, LR_diff (mean_L - mean_R)
Outcomes: mes_coc, mes_tot_miss, bit_tot_miss, bit_coc, nih_total
  → thresholded per outcome to avoid floor effects
Confounders: age, gender, lesion_side — residualised before Spearman
FDR-BH across all tests; scatter plots for strongest results.
"""
import pickle, warnings
import numpy as np
import pandas as pd
from scipy import signal as scipy_signal
from scipy.stats import spearmanr, zscore, f as f_dist
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

OUT = '/home/cemal/Desktop/Opus/vph/vph_updated'

# ── 1. LOAD DATA ─────────────────────────────────────────────
with open('/home/cemal/Desktop/Opus/stroke_csv_in_progress.pkl', 'rb') as f:
    df = pickle.load(f)

patients = df[(df.subj_type == 0) & df.inclusion_notes.isna()].copy()

# ── 2. FEATURE EXTRACTION — FIRST RUN ONLY ───────────────────
def spectral_entropy(x):
    x = x[np.isfinite(x)]
    if len(x) < 8:
        return np.nan
    _, pxx = scipy_signal.welch(x, nperseg=min(64, len(x)))
    pxx = pxx / (pxx.sum() + 1e-12)
    return float(-np.sum(pxx * np.log(pxx + 1e-12)))

def features_from_run(arr, run_idx=0):
    """Extract features from a single run (column) of a (TRs, runs) array."""
    if not isinstance(arr, np.ndarray) or arr.ndim != 2:
        return dict(mean=np.nan, var=np.nan, spec_ent=np.nan, ac1=np.nan)
    if arr.shape[1] <= run_idx:
        return dict(mean=np.nan, var=np.nan, spec_ent=np.nan, ac1=np.nan)
    x = arr[:, run_idx].astype(float)
    v = x[np.isfinite(x)]
    if len(v) < 10:
        return dict(mean=np.nan, var=np.nan, spec_ent=np.nan, ac1=np.nan)
    ac1 = float(np.corrcoef(v[:-1], v[1:])[0, 1]) if len(v) > 2 else np.nan
    return dict(
        mean    = float(np.mean(v)),
        var     = float(np.var(v, ddof=1)),
        spec_ent= spectral_entropy(v),
        ac1     = ac1,
    )

feat_rows = []
for _, row in patients.iterrows():
    fl = features_from_run(row.get('r_trtrcorr_eye_l'), run_idx=0)
    fr = features_from_run(row.get('r_trtrcorr_eye_r'), run_idx=0)
    # LR_diff: signed asymmetry (L - R) for mean of first run
    lr_diff = (fl['mean'] - fr['mean']) if (np.isfinite(fl['mean']) and np.isfinite(fr['mean'])) else np.nan
    # Bilateral mean for all other features
    def bil(a, b):
        v = [x for x in [a, b] if np.isfinite(x)]
        return float(np.mean(v)) if v else np.nan
    feat_rows.append({
        'ID':      row['ID'],
        'Session': row['Session'],
        'mean':    bil(fl['mean'],     fr['mean']),
        'var':     bil(fl['var'],      fr['var']),
        'spec_ent':bil(fl['spec_ent'], fr['spec_ent']),
        'ac1':     bil(fl['ac1'],      fr['ac1']),
        'lr_diff': lr_diff,
    })

feat_df = pd.DataFrame(feat_rows)
print(f"Feature rows: {len(feat_df)}  |  patients: {feat_df.ID.nunique()}")
print(feat_df.groupby('Session')[['mean','var','spec_ent','ac1','lr_diff']].count())

# ── 3. MERGE WITH OUTCOMES & CONFOUNDERS ─────────────────────
OUTCOMES = [
    ('mes_coc',      'Mesulam CoC'),
    ('mes_tot_miss', 'Mesulam misses'),
    ('bit_tot_miss', 'BIT misses'),
    ('bit_coc',      'BIT CoC'),
    ('nih_total',    'NIHSS total'),
]
OUT_COLS = [c for c, _ in OUTCOMES if c in patients.columns]
OUTCOMES = [(c, l) for c, l in OUTCOMES if c in patients.columns]

FEATURES = [
    ('mean',     'Mean r(t)'),
    ('var',      'Variance'),
    ('spec_ent', 'Spectral entropy'),
    ('ac1',      'AC1'),
    ('lr_diff',  'LR asymmetry'),
]
FEAT_COLS = [c for c, _ in FEATURES]

clin = patients[patients.Session.isin(['acute','followup','followup2'])][
    ['ID','Session'] + OUT_COLS + ['age','gender','lesion_side']].copy()

merged = clin.merge(feat_df, on=['ID','Session'], how='inner')
print(f"\nMerged rows: {len(merged)}, patients: {merged.ID.nunique()}")
print(merged.groupby('Session').size().to_dict())

# ── 4. PARTIAL SPEARMAN — ACUTE ONLY ────────────────────────
# Bilateral features: control for lesion side (not a mediator)
# LR asymmetry: lesion side is part of the causal pathway → do NOT control for it
CONFOUNDERS_BILAT = ['age', 'gender', 'lesion_side']
CONFOUNDERS_LR    = ['age', 'gender']
FEATURE_CONFOUNDERS = {
    'mean':     CONFOUNDERS_BILAT,
    'var':      CONFOUNDERS_BILAT,
    'spec_ent': CONFOUNDERS_BILAT,
    'ac1':      CONFOUNDERS_BILAT,
    'lr_diff':  CONFOUNDERS_LR,
}

def residualise(vec, conf_df):
    """OLS-residualise vec against columns of conf_df; return residuals."""
    mask = np.isfinite(vec) & conf_df.notna().all(axis=1).values
    out = vec.copy().astype(float)
    if mask.sum() < 10:
        return out
    X = conf_df.values[mask].astype(float)
    y = vec[mask]
    coef, *_ = np.linalg.lstsq(np.column_stack([X, np.ones(len(X))]), y, rcond=None)
    out[mask] = y - np.column_stack([X, np.ones(len(X))]) @ coef
    return out


# Literature-based thresholds (from R1_analyses/rq1_sensitivity_floor.py):
#   mes_coc      > 0.083  (Rorden & Karnath 2010)
#   bit_tot_miss >= 4     (Wilson et al. 1987 BIT manual)
#   nih_total    > 0      (floor only)
#   mes_tot_miss > 0      (floor only, no literature cutoff)
#   bit_coc      > 0      (floor only, no literature cutoff)
THRESHOLDS = {
    'mes_coc':      lambda x: x > 0,
    'mes_tot_miss': lambda x: x > 0,
    'bit_tot_miss': lambda x: x >= 4,
    'bit_coc':      lambda x: x > 0,
    'nih_total':    lambda x: x > 0,
}

results = []
sub = merged[merged.Session == 'acute'].copy()

# Pre-residualise each feature with its own confounder set
feat_resid = {}
for fc, _ in FEATURES:
    conf_df = sub[FEATURE_CONFOUNDERS[fc]]
    feat_resid[fc] = residualise(sub[fc].values.astype(float), conf_df)

# Pre-residualise outcomes for each confounder set
out_resid_bilat = {}
out_resid_lr    = {}
for out_col, _ in OUTCOMES:
    raw_out = sub[out_col].values.astype(float)
    out_resid_bilat[out_col] = residualise(raw_out, sub[CONFOUNDERS_BILAT])
    out_resid_lr[out_col]    = residualise(raw_out, sub[CONFOUNDERS_LR])

for out_col, out_lbl in OUTCOMES:
    raw_out = sub[out_col].values.astype(float)
    thresh_fn = THRESHOLDS.get(out_col, lambda x: x > 0)
    keep = np.isfinite(raw_out) & thresh_fn(raw_out)
    if keep.sum() < 15:
        print(f"  Skipping {out_col}: only {keep.sum()} above threshold")
        continue

    for fc, fl in FEATURES:
        feat_vec = feat_resid[fc]
        # Use matching residualised outcome
        out_r = out_resid_bilat[out_col] if fc != 'lr_diff' else out_resid_lr[out_col]
        valid = keep & np.isfinite(feat_vec)
        if valid.sum() < 15:
            continue
        x = feat_vec[valid]
        y = out_r[valid]
        n = int(valid.sum())

        # Spearman
        rho, p_spear = spearmanr(x, y)

        # Linear R²
        c1 = np.polyfit(x, y, 1)
        ss_res_lin = np.sum((y - np.polyval(c1, x)) ** 2)
        ss_tot     = np.sum((y - y.mean()) ** 2)
        r2_lin = 1 - ss_res_lin / ss_tot if ss_tot > 0 else np.nan

        # Quadratic R² + F-test (quadratic term vs linear)
        c2 = np.polyfit(x, y, 2)
        ss_res_quad = np.sum((y - np.polyval(c2, x)) ** 2)
        r2_quad = 1 - ss_res_quad / ss_tot if ss_tot > 0 else np.nan
        if n > 3 and ss_res_quad > 0:
            F_quad = ((ss_res_lin - ss_res_quad) / 1) / (ss_res_quad / (n - 3))
            p_quad = float(f_dist.sf(F_quad, 1, n - 3))
        else:
            F_quad, p_quad = np.nan, np.nan

        results.append(dict(
            feature=fl, feat_col=fc,
            outcome=out_lbl, out_col=out_col,
            n=n,
            rho=round(rho, 4),       p_spear=p_spear,
            r2_lin=round(r2_lin, 4), r2_quad=round(r2_quad, 4),
            p_quad=p_quad,
        ))

res_df = pd.DataFrame(results)
_, res_df['p_spear_fdr'], _, _ = multipletests(res_df['p_spear'], method='fdr_bh')
_, res_df['p_quad_fdr'],  _, _ = multipletests(res_df['p_quad'],  method='fdr_bh')
res_df['sig_spear'] = res_df['p_spear_fdr'].apply(
    lambda q: '***' if q<0.001 else '**' if q<0.01 else '*' if q<0.05 else '†' if q<0.10 else '')
res_df['sig_quad'] = res_df['p_quad_fdr'].apply(
    lambda q: '***' if q<0.001 else '**' if q<0.01 else '*' if q<0.05 else '†' if q<0.10 else '')
res_df = res_df.sort_values('p_spear_fdr').reset_index(drop=True)

# ── 5. PRINT RESULTS ─────────────────────────────────────────
print("\n" + "="*100)
print("SPEARMAN + QUADRATIC F-TEST | ACUTE ONLY | FDR-BH | mes_coc >0 | bit_tot_miss >=4 | others >0")
print("="*100)
print(f"\n  {'Feature':<20} {'Outcome':<22} {'n':>5} {'rho':>7} {'p_sp_fdr':>10} {'sig':>4}  "
      f"{'R²lin':>6} {'R²quad':>7} {'p_q_fdr':>9} {'sig_q':>5}")
print("  " + "-"*95)
for _, r in res_df.iterrows():
    flag_s = ' ◄' if r.p_spear_fdr < 0.05 else (' †' if r.p_spear_fdr < 0.10 else '  ')
    flag_q = ' ◄' if r.p_quad_fdr  < 0.05 else (' †' if r.p_quad_fdr  < 0.10 else '  ')
    print(f"  {r.feature:<20} {r.outcome:<22} {r.n:>5} {r.rho:>7.3f} {r.p_spear_fdr:>10.4f} "
          f"{r.sig_spear:>4}{flag_s}  "
          f"{r.r2_lin:>6.3f} {r.r2_quad:>7.3f} {r.p_quad_fdr:>9.4f} {r.sig_quad:>5}{flag_q}")

res_df.to_csv(f'{OUT}/vph_clean_correlations.csv', index=False)
print(f"\nSaved vph_clean_correlations.csv")

# ── 6. SCATTER PLOTS — polynomial fit helper ─────────────────
def poly2_r2(x, y):
    """Fit degree-2 polynomial, return coefficients and R²."""
    coeffs = np.polyfit(x, y, 2)
    y_hat  = np.polyval(coeffs, x)
    ss_res = np.sum((y - y_hat) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return coeffs, r2

def plot_scatter(ax, x, y, row, label_p, p_val):
    ax.scatter(x, y, s=35, alpha=0.6, color='steelblue', edgecolors='none')
    xl = np.linspace(x.min(), x.max(), 200)
    coeffs, r2 = poly2_r2(x, y)
    ax.plot(xl, np.polyval(coeffs, xl), color='firebrick', linewidth=2)
    sign = '+' if row['rho'] > 0 else ''
    ax.set_title(
        f"{row['feature']} → {row['outcome']}\n"
        f"ρ={sign}{row['rho']:.3f}, {label_p}={p_val:.4f}, R²={r2:.3f}, n={row['n']}",
        fontsize=9, fontweight='bold'
    )
    ax.set_xlabel(f"{row['feature']} (residualised)", fontsize=8)
    ax.set_ylabel(f"{row['outcome']} (residualised)", fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)

# ── scatter plots — top 9 strongest FDR-sig results ──────────
sig_df = res_df[(res_df.p_spear_fdr < 0.05) | (res_df.p_quad_fdr < 0.05)].head(9)
print(f"\nPlotting {len(sig_df)} significant results")

if len(sig_df) > 0:
    ncols = 3
    nrows = int(np.ceil(len(sig_df) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4.8, nrows * 4.2))
    axes = np.array(axes).flatten()

    acute_sub = merged[merged.Session == 'acute'].copy()
    for idx, (_, row) in enumerate(sig_df.iterrows()):
        ax = axes[idx]
        fc  = row['feat_col']
        oc  = row['out_col']
        raw_out = acute_sub[oc].values.astype(float)
        thresh_fn = THRESHOLDS.get(oc, lambda x: x > 0)
        keep = np.isfinite(raw_out) & thresh_fn(raw_out)
        conf_plot = acute_sub[FEATURE_CONFOUNDERS[fc]]
        x = residualise(acute_sub[fc].values.astype(float), conf_plot)
        y = residualise(raw_out, conf_plot)
        valid = keep & np.isfinite(x)
        plot_scatter(ax, x[valid], y[valid], row, 'FDR_sp', row['p_spear_fdr'])

    for idx in range(len(sig_df), len(axes)):
        axes[idx].set_visible(False)

    fig.suptitle(
        'Oculomotor features × clinical outcomes (residualised: age, gender, lesion side)\n'
        'First-run, partial Spearman, poly-2 fit',
        fontsize=13, fontweight='bold'
    )
    plt.tight_layout()
    plt.savefig(f'{OUT}/vph_scatter_sig.png', dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved vph_scatter_sig.png")
else:
    print("No FDR-sig results — saving top 9 by raw p instead")
    top_df = res_df.head(9)
    ncols = 3
    nrows = int(np.ceil(len(top_df) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4.8, nrows * 4.2))
    axes = np.array(axes).flatten()
    acute_sub = merged[merged.Session == 'acute'].copy()
    for idx, (_, row) in enumerate(top_df.iterrows()):
        ax = axes[idx]
        fc  = row['feat_col']
        oc  = row['out_col']
        raw_out = acute_sub[oc].values.astype(float)
        thresh_fn = THRESHOLDS.get(oc, lambda x: x > 0)
        keep = np.isfinite(raw_out) & thresh_fn(raw_out)
        conf_plot = acute_sub[FEATURE_CONFOUNDERS[fc]]
        x = residualise(acute_sub[fc].values.astype(float), conf_plot)
        y = residualise(raw_out, conf_plot)
        valid = keep & np.isfinite(x)
        plot_scatter(ax, x[valid], y[valid], row, 'p_sp', row['p_spear'])
    for idx in range(len(top_df), len(axes)):
        axes[idx].set_visible(False)
    fig.suptitle(
        'Oculomotor features × clinical outcomes — top 9 (residualised: age, gender, lesion side)\n'
        'First-run, partial Spearman, poly-2 fit',
        fontsize=13, fontweight='bold'
    )
    plt.tight_layout()
    plt.savefig(f'{OUT}/vph_scatter_top9.png', dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved vph_scatter_top9.png")

print("\nDone.")
