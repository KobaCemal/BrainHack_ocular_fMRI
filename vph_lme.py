"""
VPH 2026 — longitudinal LME
Does acute oculomotor feature (first run) predict recovery trajectory?

Model per feature × outcome (25 combinations):
  Bilateral (mean, var, spec_ent, ac1):
    outcome ~ log_days_z * feat_acute_z + age_z + gender + lesion_side + lesion_vol_z + (1|ID)
  LR asymmetry (lesion_side is mediator, not confounder):
    outcome ~ log_days_z * feat_acute_z + age_z + gender + lesion_vol_z + (1|ID)

Key term: log_days_z:feat_acute_z
  Negative β → higher acute feature predicts steeper improvement over time
FDR-BH across all 25 interaction p-values.
"""
import pickle, warnings
import numpy as np
import pandas as pd
from scipy import signal as scipy_signal
from scipy.stats import zscore
import statsmodels.formula.api as smf
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

# ── 2. FEATURE EXTRACTION — FIRST RUN ONLY (acute) ───────────
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
for _, row in patients.iterrows():
    fl = features_from_run(row.get('r_trtrcorr_eye_l'))
    fr = features_from_run(row.get('r_trtrcorr_eye_r'))
    def bil(a, b):
        v = [x for x in [a, b] if np.isfinite(x)]
        return float(np.mean(v)) if v else np.nan
    lr_diff = (fl.get('mean', np.nan) - fr.get('mean', np.nan))
    if not (np.isfinite(fl.get('mean', np.nan)) and np.isfinite(fr.get('mean', np.nan))):
        lr_diff = np.nan
    feat_rows.append({'ID': row['ID'], 'Session': row['Session'],
        'mean':     bil(fl.get('mean',np.nan),     fr.get('mean',np.nan)),
        'var':      bil(fl.get('var',np.nan),      fr.get('var',np.nan)),
        'spec_ent': bil(fl.get('spec_ent',np.nan), fr.get('spec_ent',np.nan)),
        'ac1':      bil(fl.get('ac1',np.nan),      fr.get('ac1',np.nan)),
        'lr_diff':  lr_diff})

feat_df = pd.DataFrame(feat_rows)
acute_feats = feat_df[feat_df.Session == 'acute'][
    ['ID','mean','var','spec_ent','ac1','lr_diff']].copy()
acute_feats.columns = ['ID'] + [f'{c}_acute' for c in ['mean','var','spec_ent','ac1','lr_diff']]

# ── 3. BUILD LONGITUDINAL TABLE ───────────────────────────────
OUTCOMES = [
    ('mes_coc',      'Mesulam CoC'),
    ('mes_tot_miss', 'Mesulam misses'),
    ('bit_tot_miss', 'BIT misses'),
    ('bit_coc',      'BIT CoC'),
    ('nih_total',    'NIHSS total'),
]
OUTCOMES = [(c, l) for c, l in OUTCOMES if c in patients.columns]
OUT_COLS = [c for c, _ in OUTCOMES]

SESSION_SCAN = {'acute': 'acute_scan', 'followup': 'three_m_scan', 'followup2': 'one_y_scan'}

clin = patients[patients.Session.isin(SESSION_SCAN)][
    ['ID','Session'] + OUT_COLS +
    ['age','gender','lesion_side','lesion_volume_mm3','date_stroke',
     'acute_scan','three_m_scan','one_y_scan']].copy()

# Compute actual days since stroke onset; cap followup2 at 450 (scheduling outliers)
def get_days(row):
    scan_col = SESSION_SCAN.get(row['Session'])
    if scan_col is None: return np.nan
    try:
        d = (pd.to_datetime(row[scan_col]) - pd.to_datetime(row['date_stroke'])).days
        if row['Session'] == 'followup2' and d > 450:
            d = 450
        return float(d)
    except Exception:
        return np.nan

clin['days'] = clin.apply(get_days, axis=1)
clin = clin[clin.days.notna() & (clin.days > 0)].copy()
clin['log_days'] = np.log(clin['days'])

# Propagate acute lesion volume to follow-up rows
acute_lv = (clin[clin.Session == 'acute'][['ID','lesion_volume_mm3']]
            .rename(columns={'lesion_volume_mm3': '_lv_acute'}))
clin = clin.merge(acute_lv, on='ID', how='left')
clin['lesion_volume_mm3'] = clin['lesion_volume_mm3'].fillna(clin['_lv_acute'])

# Merge acute features
long_df = clin.merge(acute_feats, on='ID', how='inner')
print(f"Longitudinal table: {long_df.ID.nunique()} subjects, {len(long_df)} rows")
print(long_df.groupby('Session').size().to_dict())

# Merge acute outcome scores as fixed baseline covariates (for baseline-controlled models)
acute_scores = (long_df[long_df.Session == 'acute'][['ID'] + OUT_COLS]
                .rename(columns={c: f'{c}_acute_score' for c in OUT_COLS}))
long_df = long_df.merge(acute_scores, on='ID', how='left')

# ── 4. Z-SCORE ────────────────────────────────────────────────
# log_days: z-score across all rows
long_df['log_days_z'] = zscore(long_df['log_days'], nan_policy='omit')

# Continuous covariates: z-score on acute values, apply same scale to all rows
for col in ['age', 'lesion_volume_mm3']:
    acute_vals = long_df.loc[long_df.Session == 'acute', col]
    mu, sd = acute_vals.mean(), acute_vals.std()
    long_df[f'{col}_z'] = (long_df[col] - mu) / sd

# Features: z-score using acute mean/SD
FEAT_COLS = ['mean_acute','var_acute','spec_ent_acute','ac1_acute','lr_diff_acute']
for fc in FEAT_COLS:
    mu, sd = long_df[fc].mean(), long_df[fc].std()
    long_df[f'{fc}_z'] = (long_df[fc] - mu) / (sd if sd > 0 else 1)

# Baseline outcome scores: z-score
for oc, _ in OUTCOMES:
    col = f'{oc}_acute_score'
    mu, sd = long_df[col].mean(), long_df[col].std()
    long_df[f'{col}_z'] = (long_df[col] - mu) / (sd if sd > 0 else 1)

# ── 5. LME MODELS ─────────────────────────────────────────────
FEATURES = [
    ('mean_acute_z',     'Mean r(t)',       False),  # (col, label, is_lr)
    ('var_acute_z',      'Variance',        False),
    ('spec_ent_acute_z', 'Spectral entropy',False),
    ('ac1_acute_z',      'AC1',             False),
    ('lr_diff_acute_z',  'LR asymmetry',    True),
]

results = []
print("\nRunning LME models...")

for feat_col, feat_lbl, is_lr in FEATURES:
    for out_col, out_lbl in OUTCOMES:
        # Build analysis subset
        covs = [out_col, feat_col, 'log_days_z', 'age_z', 'gender',
                'lesion_volume_mm3_z', 'lesion_side']
        if is_lr:
            covs = [c for c in covs if c != 'lesion_side']

        sub = long_df[['ID'] + covs].dropna().copy()
        if sub.ID.nunique() < 20 or len(sub) < 40:
            print(f"  Skip {feat_lbl} × {out_lbl}: n={sub.ID.nunique()} subjects")
            continue

        # Interaction term
        sub['feat_x_time'] = sub[feat_col] * sub['log_days_z']

        # Formula
        cov_parts = ['log_days_z', feat_col, 'feat_x_time',
                     'age_z', 'gender', 'lesion_volume_mm3_z']
        if not is_lr:
            cov_parts.append('lesion_side')
        formula = f'{out_col} ~ ' + ' + '.join(cov_parts)

        sub = sub.sort_values(['ID', 'log_days_z'])
        try:
            mdl = smf.mixedlm(formula, data=sub, groups=sub['ID'])
            res = mdl.fit(reml=True, method='lbfgs')
        except Exception:
            try:
                res = mdl.fit(reml=False, method='nm')
            except Exception as e:
                print(f"  FAIL {feat_lbl} × {out_lbl}: {e}")
                continue

        if 'feat_x_time' not in res.params:
            continue

        b_int  = res.params['feat_x_time']
        p_int  = res.pvalues['feat_x_time']
        b_main = res.params.get(feat_col, np.nan)
        p_main = res.pvalues.get(feat_col, np.nan)
        b_time = res.params.get('log_days_z', np.nan)
        p_time = res.pvalues.get('log_days_z', np.nan)

        results.append(dict(
            feature=feat_lbl, feat_col=feat_col, is_lr=is_lr,
            outcome=out_lbl,  out_col=out_col,
            n_subj=sub.ID.nunique(), n_obs=len(sub),
            b_int=round(b_int,5),   p_int=p_int,
            b_main=round(b_main,5), p_main=p_main,
            b_time=round(b_time,5), p_time=p_time,
        ))

res_df = pd.DataFrame(results)
_, res_df['p_fdr'], _, _ = multipletests(res_df['p_int'], method='fdr_bh')
res_df['sig'] = res_df['p_fdr'].apply(
    lambda q: '***' if q<0.001 else '**' if q<0.01 else '*' if q<0.05 else '†' if q<0.10 else '')
res_df = res_df.sort_values('p_fdr').reset_index(drop=True)

# ── 6. PRINT RESULTS ──────────────────────────────────────────
print("\n" + "="*90)
print("LME: outcome ~ log_days_z * feat_acute_z + age + gender + lesion_vol [+ lesion_side]")
print("KEY TERM: feat_x_time  |  FDR-BH across all 25 interaction terms")
print("="*90)
print(f"\n  {'Feature':<20} {'Outcome':<22} {'n_subj':>6} {'β_int':>9} {'p_int':>8} {'FDR':>8} {'sig'}")
print("  " + "-"*80)
for _, r in res_df.iterrows():
    flag = ' ◄' if r.p_fdr < 0.05 else (' †' if r.p_fdr < 0.10 else '')
    print(f"  {r.feature:<20} {r.outcome:<22} {r.n_subj:>6} "
          f"{r.b_int:>9.4f} {r.p_int:>8.4f} {r.p_fdr:>8.4f} {r.sig}{flag}")

res_df.to_csv(f'{OUT}/vph_lme_results.csv', index=False)
print(f"\nSaved vph_lme_results.csv")

# ── 6b. BASELINE-CONTROLLED MODELS ───────────────────────────
print("\n" + "="*90)
print("LME + BASELINE OUTCOME: same model + acute outcome score as covariate")
print("Tests whether feature predicts recovery BEYOND baseline severity")
print("="*90)

results_bc = []
for feat_col, feat_lbl, is_lr in FEATURES:
    for out_col, out_lbl in OUTCOMES:
        baseline_col   = f'{out_col}_acute_score_z'
        covs = [out_col, feat_col, 'log_days_z', 'age_z', 'gender',
                'lesion_volume_mm3_z', 'lesion_side', baseline_col]
        if is_lr:
            covs = [c for c in covs if c != 'lesion_side']

        sub = long_df[['ID'] + covs].dropna().copy()
        if sub.ID.nunique() < 20 or len(sub) < 40:
            continue

        sub['feat_x_time'] = sub[feat_col] * sub['log_days_z']
        cov_parts = ['log_days_z', feat_col, 'feat_x_time',
                     'age_z', 'gender', 'lesion_volume_mm3_z', baseline_col]
        if not is_lr:
            cov_parts.append('lesion_side')
        formula = f'{out_col} ~ ' + ' + '.join(cov_parts)

        sub = sub.sort_values(['ID', 'log_days_z'])
        try:
            mdl = smf.mixedlm(formula, data=sub, groups=sub['ID'])
            res = mdl.fit(reml=True, method='lbfgs')
        except Exception:
            try:
                res = mdl.fit(reml=False, method='nm')
            except Exception as e:
                print(f"  FAIL {feat_lbl} × {out_lbl}: {e}")
                continue

        if 'feat_x_time' not in res.params:
            continue

        b_int = res.params['feat_x_time']
        p_int = res.pvalues['feat_x_time']
        results_bc.append(dict(
            feature=feat_lbl, feat_col=feat_col, is_lr=is_lr,
            outcome=out_lbl,  out_col=out_col,
            n_subj=sub.ID.nunique(), n_obs=len(sub),
            b_int=round(b_int,5), p_int=p_int,
        ))

res_bc_df = pd.DataFrame(results_bc)
_, res_bc_df['p_fdr'], _, _ = multipletests(res_bc_df['p_int'], method='fdr_bh')
res_bc_df['sig'] = res_bc_df['p_fdr'].apply(
    lambda q: '***' if q<0.001 else '**' if q<0.01 else '*' if q<0.05 else '†' if q<0.10 else '')
res_bc_df = res_bc_df.sort_values('p_fdr').reset_index(drop=True)

print(f"\n  {'Feature':<20} {'Outcome':<22} {'n_subj':>6} {'β_int':>9} {'p_int':>8} {'FDR':>8} {'sig'}")
print("  " + "-"*80)
for _, r in res_bc_df.iterrows():
    flag = ' ◄' if r.p_fdr < 0.05 else (' †' if r.p_fdr < 0.10 else '')
    print(f"  {r.feature:<20} {r.outcome:<22} {r.n_subj:>6} "
          f"{r.b_int:>9.4f} {r.p_int:>8.4f} {r.p_fdr:>8.4f} {r.sig}{flag}")

res_bc_df.to_csv(f'{OUT}/vph_lme_results_baseline_ctrl.csv', index=False)
print(f"\nSaved vph_lme_results_baseline_ctrl.csv")

# ── 7. HEATMAP ────────────────────────────────────────────────
FEAT_ORDER = ['Mean r(t)', 'Variance', 'Spectral entropy', 'AC1', 'LR asymmetry']
OUT_ORDER  = [l for _, l in OUTCOMES]

pivot_b   = res_df.pivot(index='feature', columns='outcome', values='b_int') \
                  .reindex(index=FEAT_ORDER, columns=OUT_ORDER)
pivot_fdr = res_df.pivot(index='feature', columns='outcome', values='p_fdr') \
                  .reindex(index=FEAT_ORDER, columns=OUT_ORDER)

fig, ax = plt.subplots(figsize=(11, 4.5))
bvals = pivot_b.values.astype(float)
vlim  = np.nanmax(np.abs(bvals)) * 1.05
im = ax.imshow(bvals, cmap='RdBu_r', vmin=-vlim, vmax=vlim, aspect='auto')
ax.set_xticks(range(len(OUT_ORDER)));  ax.set_xticklabels(OUT_ORDER, rotation=30, ha='right', fontsize=10)
ax.set_yticks(range(len(FEAT_ORDER))); ax.set_yticklabels(FEAT_ORDER, fontsize=10)
for i in range(len(FEAT_ORDER)):
    for j in range(len(OUT_ORDER)):
        b   = pivot_b.values[i, j]
        fdr = pivot_fdr.values[i, j]
        if not np.isfinite(b): continue
        star = '***' if fdr<0.001 else '**' if fdr<0.01 else '*' if fdr<0.05 else '†' if fdr<0.10 else ''
        tc = 'white' if abs(b) > vlim * 0.55 else 'black'
        ax.text(j, i, f'{b:.3f}\n{star}', ha='center', va='center', fontsize=8.5, color=tc)
plt.colorbar(im, ax=ax, label='β (interaction: feat × log_days)')
ax.set_title('LME interaction β: acute oculomotor feature × log(days)\n'
             'Negative β = higher acute feature → steeper recovery  |  FDR-BH corrected',
             fontsize=11, fontweight='bold')
plt.tight_layout()
plt.savefig(f'{OUT}/vph_lme_heatmap.png', dpi=150, bbox_inches='tight', facecolor='white')
plt.close()
print("Saved vph_lme_heatmap.png")

# ── 8. TRAJECTORY PLOTS — FDR-significant interactions ────────
sig_rows = res_df[res_df.p_fdr < 0.10].head(6)   # up to 6 panels (trends included)
print(f"\nTrajectory plots: {len(sig_rows)} panels (FDR < 0.10)")

TERT_COLORS = {'Low': '#2196F3', 'Mid': '#9E9E9E', 'High': '#F44336'}
TERT_ORDER  = ['Low', 'Mid', 'High']

if len(sig_rows) > 0:
    ncols = min(3, len(sig_rows))
    nrows = int(np.ceil(len(sig_rows) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 5.5, nrows * 4.8), sharey=False)
    axes = np.array(axes).flatten()

    for ax_idx, (_, row) in enumerate(sig_rows.iterrows()):
        ax = axes[ax_idx]
        fc  = row['feat_col']
        oc  = row['out_col']
        is_lr = row['is_lr']

        # Subset with required columns
        need = ['ID','Session','days','log_days_z', fc, oc]
        sub = long_df[need].dropna().copy()

        # Tertile split on acute feature value
        acute_vals = sub[sub.Session == 'acute'].set_index('ID')[fc]
        if len(acute_vals) < 9:
            axes[ax_idx].set_visible(False)
            continue
        q33 = float(acute_vals.quantile(0.33))
        q67 = float(acute_vals.quantile(0.67))
        subj_tert = acute_vals.apply(
            lambda x: 'Low' if x <= q33 else ('High' if x >= q67 else 'Mid'))
        sub['tertile'] = sub['ID'].map(subj_tert)
        sub = sub.dropna(subset=['tertile'])
        n_by_tert = subj_tert.value_counts()

        # Group mean ± SEM
        for tert in TERT_ORDER:
            color = TERT_COLORS[tert]
            n_t   = n_by_tert.get(tert, 0)
            xs, ms, sems = [], [], []
            for ses in ['acute','followup','followup2']:
                vals = sub[(sub.Session == ses) & (sub.tertile == tert)][oc].dropna()
                if len(vals) >= 3:
                    day_med = sub[sub.Session == ses]['days'].median()
                    xs.append(day_med); ms.append(vals.mean()); sems.append(vals.sem())
            if not ms: continue
            xs = np.array(xs); ms = np.array(ms); sems = np.array(sems)
            ax.fill_between(xs, ms-sems, ms+sems, color=color, alpha=0.18, zorder=2)
            ax.plot(xs, ms, color=color, linewidth=2.5, marker='o',
                    markersize=7, zorder=3, label=f'{tert} (n={n_t})')

        ax.set_xscale('log')
        ax.set_xticks([14, 90, 380])
        ax.set_xticklabels(['Acute\n(~2 wk)', '3 months', '1 year'], fontsize=9)
        ax.set_xlabel('Days post-stroke (log scale)', fontsize=9)
        ax.set_ylabel(row['outcome'], fontsize=10)
        sign = '↓' if row['b_int'] < 0 else '↑'
        ax.set_title(f"{row['feature']} × {row['outcome']}\n"
                     f"β={row['b_int']:.3f}, FDR={row['p_fdr']:.4f} {row['sig']}",
                     fontsize=10, fontweight='bold')
        ax.legend(title=f"Acute {row['feature']} tertile",
                  fontsize=8, title_fontsize=8, framealpha=0.6)
        ax.text(0.02, 0.02,
                f"Low ≤{q33:.3f}  |  Mid  |  High ≥{q67:.3f}",
                transform=ax.transAxes, fontsize=7.5, color='#444444',
                va='bottom', ha='left')
        ax.spines[['top','right']].set_visible(False)

    for ax_idx in range(len(sig_rows), len(axes)):
        axes[ax_idx].set_visible(False)

    fig.suptitle('LME trajectory plots — acute oculomotor feature tertiles\n'
                 'Thin lines = individual  |  Thick = group mean ± SEM  |  '
                 'Blue=Low, Grey=Mid, Red=High',
                 fontsize=12, fontweight='bold', y=1.01)
    plt.tight_layout()
    plt.savefig(f'{OUT}/vph_lme_trajectories.png', dpi=150,
                bbox_inches='tight', facecolor='white')
    plt.close()
    print("Saved vph_lme_trajectories.png")
else:
    print("No FDR<0.10 results — skipping trajectory plots")

# ── 9. BASELINE-CONTROLLED HEATMAP ───────────────────────────
pivot_b_bc  = res_bc_df.pivot(index='feature', columns='outcome', values='b_int') \
                       .reindex(index=FEAT_ORDER, columns=OUT_ORDER)
pivot_fdr_bc= res_bc_df.pivot(index='feature', columns='outcome', values='p_fdr') \
                       .reindex(index=FEAT_ORDER, columns=OUT_ORDER)

fig, ax = plt.subplots(figsize=(11, 4.5))
bvals = pivot_b_bc.values.astype(float)
vlim  = np.nanmax(np.abs(bvals)) * 1.05 if np.any(np.isfinite(bvals)) else 1
im = ax.imshow(bvals, cmap='RdBu_r', vmin=-vlim, vmax=vlim, aspect='auto')
ax.set_xticks(range(len(OUT_ORDER)));  ax.set_xticklabels(OUT_ORDER, rotation=30, ha='right', fontsize=10)
ax.set_yticks(range(len(FEAT_ORDER))); ax.set_yticklabels(FEAT_ORDER, fontsize=10)
for i in range(len(FEAT_ORDER)):
    for j in range(len(OUT_ORDER)):
        b   = pivot_b_bc.values[i, j]
        fdr = pivot_fdr_bc.values[i, j]
        if not np.isfinite(b): continue
        star = '***' if fdr<0.001 else '**' if fdr<0.01 else '*' if fdr<0.05 else '†' if fdr<0.10 else ''
        tc = 'white' if abs(b) > vlim * 0.55 else 'black'
        ax.text(j, i, f'{b:.3f}\n{star}', ha='center', va='center', fontsize=8.5, color=tc)
plt.colorbar(im, ax=ax, label='β (interaction: feat × log_days)')
ax.set_title('LME interaction β — BASELINE OUTCOME CONTROLLED\n'
             'outcome ~ feat×log_days + acute_outcome + age + gender + lesion_vol [+lesion_side]\n'
             'Negative β = higher acute feature → steeper recovery  |  FDR-BH corrected',
             fontsize=10, fontweight='bold')
plt.tight_layout()
plt.savefig(f'{OUT}/vph_lme_heatmap_baseline_ctrl.png', dpi=150, bbox_inches='tight', facecolor='white')
plt.close()
print("Saved vph_lme_heatmap_baseline_ctrl.png")

# ── 10. BASELINE-CONTROLLED TRAJECTORY PLOTS ─────────────────
sig_rows_bc = res_bc_df[res_bc_df.p_fdr < 0.10].head(6)
print(f"\nBaseline-controlled trajectory plots: {len(sig_rows_bc)} panels (FDR < 0.10)")

if len(sig_rows_bc) > 0:
    ncols = min(3, len(sig_rows_bc))
    nrows = int(np.ceil(len(sig_rows_bc) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 5.5, nrows * 4.8), sharey=False)
    axes = np.array(axes).flatten()

    for ax_idx, (_, row) in enumerate(sig_rows_bc.iterrows()):
        ax = axes[ax_idx]
        fc  = row['feat_col']
        oc  = row['out_col']

        need = ['ID','Session','days', fc, oc]
        sub = long_df[need].dropna().copy()

        acute_vals = sub[sub.Session == 'acute'].set_index('ID')[fc]
        if len(acute_vals) < 9:
            axes[ax_idx].set_visible(False)
            continue
        q33 = float(acute_vals.quantile(0.33))
        q67 = float(acute_vals.quantile(0.67))
        subj_tert = acute_vals.apply(
            lambda x: 'Low' if x <= q33 else ('High' if x >= q67 else 'Mid'))
        sub['tertile'] = sub['ID'].map(subj_tert)
        sub = sub.dropna(subset=['tertile'])
        n_by_tert = subj_tert.value_counts()

        for tert in TERT_ORDER:
            color = TERT_COLORS[tert]
            n_t   = n_by_tert.get(tert, 0)
            xs, ms, sems = [], [], []
            for ses in ['acute','followup','followup2']:
                vals = sub[(sub.Session == ses) & (sub.tertile == tert)][oc].dropna()
                if len(vals) >= 3:
                    day_med = sub[sub.Session == ses]['days'].median()
                    xs.append(day_med); ms.append(vals.mean()); sems.append(vals.sem())
            if not ms: continue
            xs = np.array(xs); ms = np.array(ms); sems = np.array(sems)
            ax.fill_between(xs, ms-sems, ms+sems, color=color, alpha=0.18, zorder=2)
            ax.plot(xs, ms, color=color, linewidth=2.5, marker='o',
                    markersize=7, zorder=3, label=f'{tert} (n={n_t})')

        ax.set_xscale('log')
        ax.set_xticks([14, 90, 380])
        ax.set_xticklabels(['Acute\n(~2 wk)', '3 months', '1 year'], fontsize=9)
        ax.set_xlabel('Days post-stroke (log scale)', fontsize=9)
        ax.set_ylabel(row['outcome'], fontsize=10)
        ax.set_title(f"{row['feature']} × {row['outcome']}\n"
                     f"β={row['b_int']:.3f}, FDR={row['p_fdr']:.4f} {row['sig']}\n"
                     f"[baseline {row['outcome']} controlled]",
                     fontsize=9.5, fontweight='bold')
        ax.legend(title=f"Acute {row['feature']} tertile",
                  fontsize=8, title_fontsize=8, framealpha=0.6)
        ax.text(0.02, 0.02,
                f"Low ≤{q33:.3f}  |  Mid  |  High ≥{q67:.3f}",
                transform=ax.transAxes, fontsize=7.5, color='#444444',
                va='bottom', ha='left')
        ax.spines[['top','right']].set_visible(False)

    for ax_idx in range(len(sig_rows_bc), len(axes)):
        axes[ax_idx].set_visible(False)

    fig.suptitle('LME trajectories — BASELINE OUTCOME CONTROLLED\n'
                 'Interaction survives after accounting for acute severity  |  '
                 'Blue=Low, Grey=Mid, Red=High',
                 fontsize=12, fontweight='bold', y=1.01)
    plt.tight_layout()
    plt.savefig(f'{OUT}/vph_lme_trajectories_baseline_ctrl.png', dpi=150,
                bbox_inches='tight', facecolor='white')
    plt.close()
    print("Saved vph_lme_trajectories_baseline_ctrl.png")
else:
    print("No FDR<0.10 results after baseline control — skipping trajectory plots")

print("\nDone.")
