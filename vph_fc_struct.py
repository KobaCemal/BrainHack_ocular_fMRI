"""
VPH 2026 — Step 4: Oculomotor metrics × structural & functional connectivity
Part A: 5 metrics × 8 structural measures (lesion vol, disconnectome, ADC, FA)
Part B: 5 metrics × 28 FC network pairs — acute, cross-sectional
        Run for both corrmats (uncorrected) and corrmats_corrected — first run only
        FD included as covariate (new vs R4)
Part C: 5 metrics × 7 within-network FC — longitudinal LME
        FD as time-varying covariate
Method: partial Spearman (rank-residual). Confounders: age, gender, lesion_side
        (lesion_side excluded for lr_diff — causal mediator, not confounder).
        FC analyses additionally control for mean_FD and lesion_volume_mm3.
FDR-BH within each part separately.
"""
import pickle, warnings, re
import numpy as np
import pandas as pd
from scipy import stats, signal as scipy_signal
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf
from sklearn.linear_model import LinearRegression
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

OUT      = '/home/cemal/Desktop/Opus/vph'
MAP_PATH = '/home/cemal/Desktop/Opus/R4_analyses/schaefer400_network_mapping.csv'
CAP_DAYS = 450

NETWORKS = ['Vis', 'SomMot', 'DorsAttn', 'SalVentAttn', 'Limbic', 'Cont', 'Default']
NET_LABELS = {
    'Vis': 'Visual', 'SomMot': 'Somatomotor', 'DorsAttn': 'Dorsal Attn',
    'SalVentAttn': 'Ventral Attn', 'Limbic': 'Limbic',
    'Cont': 'Control', 'Default': 'Default Mode',
}
NET_COLORS = {
    'Vis': '#9467bd', 'SomMot': '#2ca02c', 'DorsAttn': '#17becf',
    'SalVentAttn': '#e377c2', 'Limbic': '#bcbd22',
    'Cont': '#ff7f0e', 'Default': '#d62728',
}

STRUCT_COLS = [
    'lesion_volume_mm3', 'disconnectome_nvox_thr_0p5',
    'fa_mean_wholebrain', 'adc_mean_wholebrain',
    'fa_mean_lesion',     'adc_mean_lesion',
    'fa_mean_disconnectome', 'adc_mean_disconnectome',
]
STRUCT_LABELS = {
    'lesion_volume_mm3':          'Lesion volume (mm³)',
    'disconnectome_nvox_thr_0p5': 'Disconnectome (vox)',
    'fa_mean_wholebrain':         'FA — whole brain',
    'adc_mean_wholebrain':        'ADC — whole brain',
    'fa_mean_lesion':             'FA — lesion',
    'adc_mean_lesion':            'ADC — lesion',
    'fa_mean_disconnectome':      'FA — disconnectome',
    'adc_mean_disconnectome':     'ADC — disconnectome',
}

METRICS = ['mean', 'var', 'spec_ent', 'ac1', 'lr_diff']
METRIC_LABELS = {
    'mean': 'Mean r(t)', 'var': 'Variance',
    'spec_ent': 'Spectral entropy', 'ac1': 'AC1', 'lr_diff': 'LR asymmetry',
}

# ── 1. LOAD ────────────────────────────────────────────────────────────────────
with open('/home/cemal/Desktop/Opus/stroke_csv_in_progress.pkl', 'rb') as f:
    df = pickle.load(f)
patients = df[(df.subj_type == 0) & df.inclusion_notes.isna()].copy()
print(f'Patients: {patients.ID.nunique()} subjects, {len(patients)} session-rows')

# ── 2. FEATURE EXTRACTION (first run, bilateral mean) ─────────────────────────
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
    lr_diff = fl.get('mean', np.nan) - fr.get('mean', np.nan)
    if not (np.isfinite(fl.get('mean', np.nan)) and np.isfinite(fr.get('mean', np.nan))):
        lr_diff = np.nan
    feat_rows.append({'ID': row['ID'], 'Session': row['Session'],
        'mean':     bil(fl.get('mean',np.nan),     fr.get('mean',np.nan)),
        'var':      bil(fl.get('var',np.nan),       fr.get('var',np.nan)),
        'spec_ent': bil(fl.get('spec_ent',np.nan),  fr.get('spec_ent',np.nan)),
        'ac1':      bil(fl.get('ac1',np.nan),       fr.get('ac1',np.nan)),
        'lr_diff':  lr_diff})

feat_df = pd.DataFrame(feat_rows)
patients = patients.merge(feat_df, on=['ID','Session'], how='left')

# ── 3. SCALAR FD (nanmean of per-run list) ─────────────────────────────────────
def scalar_fd(x):
    try:
        arr = np.array(x, dtype=float)
        return float(np.nanmean(arr))
    except Exception:
        return np.nan

patients['fd'] = patients['mean_FD'].apply(scalar_fd)

# ── 4. DAYS POST ONSET ────────────────────────────────────────────────────────
for col in ['acute_scan','three_m_scan','one_y_scan','date_stroke']:
    patients[col] = pd.to_datetime(patients[col], errors='coerce')

SESSION_SCAN = {'acute': 'acute_scan', 'followup': 'three_m_scan', 'followup2': 'one_y_scan'}
def get_scan_date(row):
    col = SESSION_SCAN.get(row['Session'])
    return row[col] if col else pd.NaT

patients['scan_date']       = patients.apply(get_scan_date, axis=1)
patients['days_post_onset'] = (patients['scan_date'] - patients['date_stroke']).dt.days
patients.loc[patients['days_post_onset'] > CAP_DAYS, 'days_post_onset'] = CAP_DAYS
patients['log_days']        = np.log(patients['days_post_onset'].clip(lower=1))

# ── 5. NETWORK FC EXTRACTION ──────────────────────────────────────────────────
print('Extracting network FC...')
net_map = pd.read_csv(MAP_PATH)
net_idx = {net: net_map[net_map['network'] == net]['roi_idx'].values - 1
           for net in NETWORKS}

FC_PAIRS  = [f'fc_{n1}_{n2}' for i, n1 in enumerate(NETWORKS)
             for j, n2 in enumerate(NETWORKS) if j >= i]
WITHIN_FC = [f'fc_{n}_{n}' for n in NETWORKS]

def extract_network_fc(mat_raw):
    if mat_raw is None or (isinstance(mat_raw, float) and np.isnan(mat_raw)):
        return {p: np.nan for p in FC_PAIRS}
    arr = np.array(mat_raw, dtype=float)
    # 3D (400,400,n_runs) → take first run; 2D already averaged
    if arr.ndim == 3:
        mat = arr[:, :, 0]
    elif arr.ndim == 2:
        mat = arr
    else:
        return {p: np.nan for p in FC_PAIRS}
    if mat.shape != (400, 400):
        return {p: np.nan for p in FC_PAIRS}
    np.fill_diagonal(mat, np.nan)
    out = {}
    for i, n1 in enumerate(NETWORKS):
        for j, n2 in enumerate(NETWORKS):
            if j < i: continue
            idx1, idx2 = net_idx[n1], net_idx[n2]
            if i == j:
                sub = mat[np.ix_(idx1, idx2)]
                out[f'fc_{n1}_{n2}'] = np.nanmean(sub[np.triu_indices(len(idx1), k=1)])
            else:
                out[f'fc_{n1}_{n2}'] = np.nanmean(mat[np.ix_(idx1, idx2)])
    return out

for fc_ver, col in [('unc', 'corrmats'), ('cor', 'corrmats_corrected')]:
    fc_rows = []
    for _, row in patients.iterrows():
        fc = extract_network_fc(row.get(col))
        fc['ID']      = row['ID']
        fc['Session'] = row['Session']
        fc_rows.append(fc)
    fc_df = pd.DataFrame(fc_rows).rename(columns={p: f'{fc_ver}_{p}' for p in FC_PAIRS})
    patients = patients.merge(fc_df, on=['ID','Session'], how='left')

print(f'FC done. Patients shape: {patients.shape}')

acute = patients[patients.Session == 'acute'].copy()
for col in ['age','gender','lesion_side','lesion_volume_mm3','fd'] + STRUCT_COLS + METRICS:
    acute[col] = pd.to_numeric(acute[col], errors='coerce')

# ── 6. UTILITIES ──────────────────────────────────────────────────────────────
def sig_stars(p):
    if pd.isna(p): return ''
    return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else '†' if p < 0.10 else ''

def partial_spearman(x_vals, y_vals, cov_list):
    """Rank-residual partial Spearman. cov_list = list of 1-D arrays."""
    cols = {'x': x_vals, 'y': y_vals}
    cols.update({f'c{i}': c for i, c in enumerate(cov_list)})
    d = pd.DataFrame(cols).replace([np.inf, -np.inf], np.nan).dropna()
    if len(d) < 10:
        return np.nan, np.nan, len(d)
    C = d[[f'c{i}' for i in range(len(cov_list))]].values
    def rank_resid(v):
        rv = stats.rankdata(v).astype(float)
        return rv - LinearRegression().fit(C, rv).predict(C)
    r, p = stats.pearsonr(rank_resid(d['x'].values), rank_resid(d['y'].values))
    return float(r), float(p), len(d)

# Confounder arrays for each metric (struct / cross-sectional)
CONF_STRUCT_BASE = ['age', 'gender']
def struct_covs(metric):
    cols = CONF_STRUCT_BASE + (['lesion_side'] if metric != 'lr_diff' else [])
    return [acute[c].values for c in cols]

# Confounder arrays for FC (add lesion_vol and FD)
CONF_FC_BASE = ['age', 'gender', 'lesion_volume_mm3', 'fd']
def fc_covs(metric):
    cols = CONF_FC_BASE + (['lesion_side'] if metric != 'lr_diff' else [])
    return [acute[c].values for c in cols]


# ══════════════════════════════════════════════════════════════════════════════
# PART A — STRUCTURAL CORRELATIONS
# ══════════════════════════════════════════════════════════════════════════════
print('\n' + '='*70)
print('PART A: STRUCTURAL CORRELATIONS (acute, partial Spearman)')
print('='*70)

struct_records = []
for metric in METRICS:
    covs = struct_covs(metric)
    for scol in STRUCT_COLS:
        r, p, n = partial_spearman(
            acute[metric].values,
            pd.to_numeric(acute[scol], errors='coerce').values,
            covs)
        struct_records.append({
            'metric': metric, 'struct': scol,
            'metric_lbl': METRIC_LABELS[metric],
            'struct_lbl': STRUCT_LABELS[scol],
            'r': r, 'p': p, 'n': n,
            'primary': scol in ['lesion_volume_mm3','disconnectome_nvox_thr_0p5'],
        })

struct_df = pd.DataFrame(struct_records)
valid_s = struct_df['p'].notna()
_, p_fdr_s, _, _ = multipletests(struct_df.loc[valid_s,'p'], method='fdr_bh')
struct_df.loc[valid_s,'p_fdr'] = p_fdr_s
struct_df['stars'] = struct_df['p_fdr'].apply(sig_stars)

print(f"\n{'Metric':<18} {'Structural':<28} {'r':>7} {'p_raw':>8} {'p_FDR':>8} {'N':>4}")
print('-'*70)
for _, row in struct_df.sort_values('p_fdr').iterrows():
    if pd.isna(row['r']): continue
    flag = ' ◄' if row.get('p_fdr',1) < 0.05 else ''
    print(f"{row['metric_lbl']:<18} {row['struct_lbl']:<28} "
          f"{row['r']:>7.3f} {row['p']:>8.4f} {row.get('p_fdr',np.nan):>8.4f} "
          f"{row['stars']:>4}{flag}  n={int(row['n'])}")

struct_df.to_csv(f'{OUT}/vph_struct_correlations.csv', index=False, float_format='%.4f')
print('\nSaved vph_struct_correlations.csv')

# ── Figure A: heatmap (5 metrics × 8 struct) ──────────────────────────────────
r_mat_s = np.full((len(METRICS), len(STRUCT_COLS)), np.nan)
p_mat_s = np.full((len(METRICS), len(STRUCT_COLS)), np.nan)
n_mat_s = np.full((len(METRICS), len(STRUCT_COLS)), np.nan)
for i, m in enumerate(METRICS):
    for j, s in enumerate(STRUCT_COLS):
        row = struct_df[(struct_df.metric==m) & (struct_df.struct==s)]
        if len(row):
            r_mat_s[i,j] = row.iloc[0]['r']
            p_mat_s[i,j] = row.iloc[0].get('p_fdr', np.nan)
            n_mat_s[i,j] = row.iloc[0]['n']

fig_a, ax_a = plt.subplots(figsize=(13, 4), constrained_layout=True)
vmax_s = max(0.4, np.nanmax(np.abs(r_mat_s)))
im_a = ax_a.imshow(r_mat_s, cmap='RdBu_r', vmin=-vmax_s, vmax=vmax_s, aspect='auto')
plt.colorbar(im_a, ax=ax_a, label='Partial Spearman r', fraction=0.02)
ax_a.set_xticks(range(len(STRUCT_COLS)))
ax_a.set_yticks(range(len(METRICS)))
ax_a.set_xticklabels([STRUCT_LABELS[s] for s in STRUCT_COLS], fontsize=8, rotation=30, ha='right')
ax_a.set_yticklabels([METRIC_LABELS[m] for m in METRICS], fontsize=9)
ax_a.axvline(1.5, color='white', lw=2.5, zorder=5)
ax_a.axvline(1.5, color='#444', lw=1.2, ls='--', zorder=6)
ax_a.text(0.5, -1.05, 'Primary', ha='center', fontsize=8, style='italic', color='#333',
          transform=ax_a.transData)
ax_a.text(4.5, -1.05, 'Diffusion (N≈28)', ha='center', fontsize=8, style='italic', color='#333',
          transform=ax_a.transData)

for i in range(len(METRICS)):
    for j in range(len(STRUCT_COLS)):
        if np.isnan(r_mat_s[i,j]): continue
        star  = sig_stars(p_mat_s[i,j])
        n_val = int(n_mat_s[i,j])
        color = 'white' if abs(r_mat_s[i,j]) > vmax_s * 0.55 else 'black'
        ax_a.text(j, i, f'{r_mat_s[i,j]:.2f}{star}\nn={n_val}',
                  ha='center', va='center', fontsize=7.5, color=color,
                  fontweight='bold' if star else 'normal')

ax_a.set_title('Oculomotor metrics × structural measures (acute)\n'
               'Partial Spearman r  |  FDR-BH (40 tests)  |  '
               'Covariates: age, gender, lesion_side (not for LR asymmetry)',
               fontsize=10, fontweight='bold')
ax_a.spines[['top','right','bottom','left']].set_visible(False)
fig_a.savefig(f'{OUT}/vph_struct_heatmap.png', dpi=200, bbox_inches='tight', facecolor='white')
plt.close()
print('Saved vph_struct_heatmap.png')

# ── Figure A2: Scatter plots — top structural correlations ─────────────────────
top_struct = struct_df.dropna(subset=['r']).sort_values('p').head(6).reset_index(drop=True)

def quad_ftest(x, y):
    """F-test for improvement of quadratic over linear fit. Returns F, p."""
    n = len(x)
    if n < 4: return np.nan, np.nan
    # Linear residuals
    c1 = np.polyfit(x, y, 1)
    ss_lin  = np.sum((y - np.polyval(c1, x))**2)
    # Quadratic residuals
    c2 = np.polyfit(x, y, 2)
    ss_quad = np.sum((y - np.polyval(c2, x))**2)
    if ss_quad == 0: return np.nan, np.nan
    F = (ss_lin - ss_quad) / 1 / (ss_quad / (n - 3))
    p = 1 - stats.f.cdf(F, 1, n - 3)
    return float(F), float(p)

ncols = 3
nrows = 2
fig_scat, axes_scat = plt.subplots(nrows, ncols, figsize=(13, 8), constrained_layout=True)
fig_scat.suptitle(
    'Oculomotor metrics × structural measures (acute)  |  Top 6 by p_raw\n'
    'Partial Spearman r (controlling for age, gender, lesion_side)\n'
    'Solid = quadratic fit  |  Dashed = linear  |  Band = 95% bootstrap CI (quadratic)',
    fontsize=10, fontweight='bold')

rng = np.random.default_rng(42)
for idx, row in top_struct.iterrows():
    ax = axes_scat.flat[idx]
    metric = row['metric']
    scol   = row['struct']

    d  = acute[[metric, scol]].apply(pd.to_numeric, errors='coerce').dropna()
    x  = d[metric].values
    y  = d[scol].values
    xs = np.linspace(x.min(), x.max(), 200)

    color = '#E64B35' if row['r'] > 0 else '#4DBBD5'
    ax.scatter(x, y, color='#444444', alpha=0.45, s=22, edgecolors='none', zorder=2)

    # Linear fit (dashed reference)
    c1 = np.polyfit(x, y, 1)
    ax.plot(xs, np.polyval(c1, xs), color=color, lw=1.5, ls='--', alpha=0.6, zorder=3)

    # Quadratic fit (solid)
    c2 = np.polyfit(x, y, 2)
    ax.plot(xs, np.polyval(c2, xs), color=color, lw=2.5, ls='-', zorder=4)

    # Bootstrap 95% CI for quadratic
    boot = []
    for _ in range(500):
        idx_b = rng.integers(0, len(x), len(x))
        cb = np.polyfit(x[idx_b], y[idx_b], 2)
        boot.append(np.polyval(cb, xs))
    boot = np.array(boot)
    ax.fill_between(xs, np.percentile(boot, 2.5, axis=0),
                    np.percentile(boot, 97.5, axis=0), color=color, alpha=0.15, zorder=1)

    # F-test for quadratic improvement
    F_val, p_quad = quad_ftest(x, y)
    quad_str = f'F-quad p={p_quad:.3f}' + (' *' if p_quad < 0.05 else (' †' if p_quad < 0.10 else ''))

    p_fdr_val = row.get('p_fdr', np.nan)
    star_str  = row['stars'] if row['stars'] else 'n.s.'
    ax.set_xlabel(row['metric_lbl'], fontsize=9)
    ax.set_ylabel(row['struct_lbl'], fontsize=9)
    ax.set_title(
        f"{row['metric_lbl']} × {row['struct_lbl']}\n"
        f"r={row['r']:+.3f}  p_FDR={p_fdr_val:.4f}{star_str}  N={int(row['n'])}  |  {quad_str}",
        fontsize=8.5, fontweight='bold' if row['p'] < 0.05 else 'normal')
    ax.tick_params(labelsize=8)
    ax.spines[['top', 'right']].set_visible(False)

fig_scat.savefig(f'{OUT}/vph_struct_scatters.png', dpi=200, bbox_inches='tight', facecolor='white')
plt.close()
print('Saved vph_struct_scatters.png')


# ══════════════════════════════════════════════════════════════════════════════
# PART B — FC CROSS-SECTIONAL (both uncorrected + lag-corrected)
# ══════════════════════════════════════════════════════════════════════════════
print('\n' + '='*70)
print('PART B: FC CORRELATIONS (acute, cross-sectional, partial Spearman)')
print('='*70)

all_fc_dfs = {}
for fc_ver, fc_name in [('unc', 'Uncorrected FC'), ('cor', 'Lag-corrected FC')]:
    print(f'\n── {fc_name} ──')
    fc_records = []
    for metric in METRICS:
        covs = fc_covs(metric)
        for pair in FC_PAIRS:
            col = f'{fc_ver}_{pair}'
            r, p, n = partial_spearman(
                acute[metric].values,
                pd.to_numeric(acute[col], errors='coerce').values if col in acute.columns else np.full(len(acute), np.nan),
                covs)
            parts = pair.replace('fc_','').split('_', 1)
            n1, n2 = parts[0], parts[1]
            fc_records.append({
                'metric': metric, 'pair': pair, 'n1': n1, 'n2': n2,
                'within': n1 == n2, 'fc_ver': fc_ver,
                'r': r, 'p': p, 'n': n,
            })

    fc_df = pd.DataFrame(fc_records)
    valid_fc = fc_df['p'].notna()
    _, p_fdr_fc, _, _ = multipletests(fc_df.loc[valid_fc,'p'], method='fdr_bh')
    fc_df.loc[valid_fc, 'p_fdr'] = p_fdr_fc
    fc_df['stars'] = fc_df['p_fdr'].apply(sig_stars)
    all_fc_dfs[fc_ver] = fc_df

    n_sig   = (fc_df['p_fdr'] < 0.05).sum()
    n_trend = ((fc_df['p_fdr'] >= 0.05) & (fc_df['p_fdr'] < 0.10)).sum()
    print(f'FDR<0.05: {n_sig}  |  Trends (0.05–0.10): {n_trend}  |  Total: {len(fc_df)} tests')
    top = fc_df.sort_values('p_fdr').head(10)
    print(f"{'Metric':<14} {'Pair':<30} {'r':>7} {'p_raw':>8} {'p_FDR':>8} {'N':>4}")
    for _, row in top.iterrows():
        if pd.isna(row['r']): continue
        flag = ' ◄' if row.get('p_fdr',1) < 0.05 else ''
        print(f"  {row['metric']:<12} {row['pair']:<30} {row['r']:>7.3f} "
              f"{row['p']:>8.4f} {row.get('p_fdr',np.nan):>8.4f} "
              f"{row['stars']:>4}{flag}  n={int(row['n'])}")

    fc_df.to_csv(f'{OUT}/vph_fc_{fc_ver}_correlations.csv', index=False, float_format='%.4f')
    print(f'Saved vph_fc_{fc_ver}_correlations.csv')

    # ── Figure B: 5-panel 7×7 heatmap ──────────────────────────────────────────
    fig_b, axes_b = plt.subplots(1, len(METRICS), figsize=(5.5*len(METRICS), 5.5),
                                  constrained_layout=True)
    fig_b.suptitle(
        f'Oculomotor metrics × {fc_name} (acute)  |  Partial Spearman r  |  FDR-BH (140 tests)\n'
        'Covariates: age, gender, lesion_vol, mean FD, lesion_side (not LR asymmetry)\n'
        'Only FDR-significant/trending cells annotated',
        fontsize=10, fontweight='bold')

    for ax_b, metric in zip(axes_b, METRICS):
        r7 = np.full((7,7), np.nan)
        p7 = np.full((7,7), np.nan)
        for i, n1 in enumerate(NETWORKS):
            for j, n2 in enumerate(NETWORKS):
                pair = f'fc_{n1}_{n2}' if j >= i else f'fc_{n2}_{n1}'
                row  = fc_df[(fc_df.metric==metric) & (fc_df.pair==pair)]
                if len(row):
                    r7[i,j] = row.iloc[0]['r']
                    r7[j,i] = row.iloc[0]['r']
                    p7[i,j] = row.iloc[0].get('p_fdr', np.nan)
                    p7[j,i] = row.iloc[0].get('p_fdr', np.nan)

        vmax_fc = max(0.25, np.nanmax(np.abs(r7)))
        im_b = ax_b.imshow(r7, cmap='RdBu_r', vmin=-vmax_fc, vmax=vmax_fc, aspect='auto')
        plt.colorbar(im_b, ax=ax_b, fraction=0.046, label='r')
        ax_b.set_xticks(range(7)); ax_b.set_yticks(range(7))
        ax_b.set_xticklabels([NET_LABELS[n] for n in NETWORKS], fontsize=6.5, rotation=40, ha='right')
        ax_b.set_yticklabels([NET_LABELS[n] for n in NETWORKS], fontsize=6.5)
        ax_b.set_title(METRIC_LABELS[metric], fontsize=10, fontweight='bold')

        # Annotate only significant/trending cells
        for i in range(7):
            for j in range(7):
                if np.isnan(r7[i,j]): continue
                star = sig_stars(p7[i,j])
                if not star: continue
                color = 'white' if abs(r7[i,j]) > vmax_fc*0.6 else 'black'
                ax_b.text(j, i, f'{r7[i,j]:.2f}{star}', ha='center', va='center',
                          fontsize=7, color=color, fontweight='bold')

        # Network colour chips
        for k, net in enumerate(NETWORKS):
            ax_b.add_patch(plt.Rectangle((-0.5, k-0.5), 0.08, 1, color=NET_COLORS[net],
                           zorder=3, transform=ax_b.get_yaxis_transform(), clip_on=False))
        ax_b.spines[['top','right','bottom','left']].set_visible(False)

    fig_b.savefig(f'{OUT}/vph_fc_{fc_ver}_heatmap.png', dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f'Saved vph_fc_{fc_ver}_heatmap.png')


# ══════════════════════════════════════════════════════════════════════════════
# PART C — LONGITUDINAL FC LME (within-network × 5 metrics)
# ══════════════════════════════════════════════════════════════════════════════
print('\n' + '='*70)
print('PART C: LONGITUDINAL FC LME (lag-corrected, within-network)')
print('='*70)

# Build long format
long = patients[patients.Session.isin(['acute','followup','followup2'])].copy()
for col in ['age','gender','lesion_side','lesion_volume_mm3','fd']:
    long[col] = pd.to_numeric(long[col], errors='coerce')

# Propagate acute metrics and lesion volume from acute session
acute_base = (patients[patients.Session=='acute']
              [['ID'] + METRICS + ['lesion_volume_mm3']]
              .rename(columns={m: f'{m}_acute' for m in METRICS})
              .rename(columns={'lesion_volume_mm3': 'lesion_vol_acute'}))
long = long.merge(acute_base, on='ID', how='inner')
print(f'Long format: {len(long)} rows, {long.ID.nunique()} subjects')

# Z-score predictors (scaled from acute-session distribution)
def z_from_acute(col_src):
    acute_vals = pd.to_numeric(long.loc[long.Session=='acute', col_src], errors='coerce').dropna()
    m_, s_ = acute_vals.mean(), acute_vals.std()
    return (pd.to_numeric(long[col_src], errors='coerce') - m_) / s_

for m in METRICS:
    long[f'{m}_z'] = z_from_acute(f'{m}_acute')
long['lesion_vol_z'] = z_from_acute('lesion_vol_acute')
long['age_z']        = z_from_acute('age')
long['log_days_z']   = z_from_acute('log_days')

# FD is time-varying — z-score from all sessions
fd_m, fd_s = long['fd'].mean(), long['fd'].std()
long['fd_z'] = (long['fd'] - fd_m) / fd_s

def _fvars(formula, data):
    tokens = re.findall(r'[A-Za-z_][A-Za-z0-9_]*', formula)
    return [t for t in tokens if t in data.columns]

def run_lme(formula, data):
    try:
        d = data.dropna(subset=_fvars(formula, data)).copy().reset_index(drop=True)
        if d.ID.nunique() < 10: return None
        md = smf.mixedlm(formula, data=d, groups=d['ID'])
        for reml in [True, False]:
            for method in ['lbfgs','bfgs','nm','powell']:
                try:
                    return md.fit(reml=reml, method=method, maxiter=500)
                except Exception:
                    continue
        return None
    except Exception as e:
        return None

lme_records = []
print(f"\n{'Metric':<18} {'Network':<16} {'β_int':>8} {'p_raw':>8}")
print('-'*55)

for metric in METRICS:
    ls_term  = '' if metric == 'lr_diff' else ' + lesion_side'
    base_cov = f'age_z + gender + lesion_vol_z{ls_term} + fd_z'

    for net in NETWORKS:
        fc_col  = f'cor_fc_{net}_{net}'
        formula = f'{fc_col} ~ log_days_z * {metric}_z + {base_cov}'
        res     = run_lme(formula, long)

        if res is None:
            lme_records.append({'metric': metric, 'network': net,
                                 'coef': np.nan, 'ci_lo': np.nan, 'ci_hi': np.nan,
                                 'p': np.nan, 'n': 0})
            continue

        # statsmodels names interaction as 'A:B' (dependent on order)
        int_term = f'log_days_z:{metric}_z'
        if int_term not in res.params.index:
            int_term = f'{metric}_z:log_days_z'
        if int_term not in res.params.index:
            int_term = next((t for t in res.params.index if metric in t and 'log' in t), None)

        if int_term is None:
            lme_records.append({'metric': metric, 'network': net,
                                 'coef': np.nan, 'ci_lo': np.nan, 'ci_hi': np.nan,
                                 'p': np.nan, 'n': int(res.nobs)})
            continue

        coef = float(res.params[int_term])
        p    = float(res.pvalues[int_term])
        try:
            ci   = res.conf_int().loc[int_term]
            ci_lo, ci_hi = float(ci.iloc[0]), float(ci.iloc[1])
        except Exception:
            ci_lo = ci_hi = np.nan

        lme_records.append({'metric': metric, 'network': net,
                             'metric_lbl': METRIC_LABELS[metric],
                             'network_lbl': NET_LABELS[net],
                             'coef': coef, 'ci_lo': ci_lo, 'ci_hi': ci_hi,
                             'p': p, 'n': int(res.nobs)})
        print(f"  {METRIC_LABELS[metric]:<16} {NET_LABELS[net]:<16} {coef:>8.4f} {p:>8.4f}{sig_stars(p)}")

lme_df = pd.DataFrame(lme_records)
valid_l = lme_df['p'].notna()
_, p_fdr_l, _, _ = multipletests(lme_df.loc[valid_l,'p'], method='fdr_bh')
lme_df.loc[valid_l,'p_fdr'] = p_fdr_l
lme_df['stars'] = lme_df['p_fdr'].apply(sig_stars)

print(f"\nFDR<0.05: {(lme_df['p_fdr']<0.05).sum()}  |  "
      f"Trends: {((lme_df['p_fdr']>=0.05)&(lme_df['p_fdr']<0.10)).sum()}  |  "
      f"Total: {valid_l.sum()} tests")

lme_df.to_csv(f'{OUT}/vph_fc_lme_results.csv', index=False, float_format='%.4f')
print('Saved vph_fc_lme_results.csv')

# ── Figure C1: LME interaction β heatmap ──────────────────────────────────────
r_lme = np.full((len(METRICS), 7), np.nan)
p_lme = np.full((len(METRICS), 7), np.nan)
n_lme = np.full((len(METRICS), 7), np.nan)

for i, m in enumerate(METRICS):
    for j, net in enumerate(NETWORKS):
        row = lme_df[(lme_df.metric==m) & (lme_df.network==net)]
        if len(row):
            r_lme[i,j] = row.iloc[0]['coef']
            p_lme[i,j] = row.iloc[0].get('p_fdr', np.nan)
            n_lme[i,j] = row.iloc[0]['n']

vmax_l = max(0.05, np.nanmax(np.abs(r_lme)))
fig_c, ax_c = plt.subplots(figsize=(11, 4), constrained_layout=True)
im_c = ax_c.imshow(r_lme, cmap='RdBu_r', vmin=-vmax_l, vmax=vmax_l, aspect='auto')
plt.colorbar(im_c, ax=ax_c, label='β  (metric × log_days)', fraction=0.02)
ax_c.set_xticks(range(7)); ax_c.set_yticks(range(len(METRICS)))
ax_c.set_xticklabels([NET_LABELS[n] for n in NETWORKS], fontsize=9, rotation=30, ha='right')
ax_c.set_yticklabels([METRIC_LABELS[m] for m in METRICS], fontsize=9)

for i in range(len(METRICS)):
    for j in range(7):
        v = r_lme[i,j]
        if np.isnan(v): continue
        star  = sig_stars(p_lme[i,j])
        n_val = int(n_lme[i,j]) if not np.isnan(n_lme[i,j]) else ''
        color = 'white' if abs(v) > vmax_l*0.6 else 'black'
        ax_c.text(j, i, f'{v:.3f}{star}\nn={n_val}',
                  ha='center', va='center', fontsize=8, color=color,
                  fontweight='bold' if star else 'normal')

ax_c.set_title(
    'Longitudinal FC LME — β for (oculomotor metric × log_days) interaction\n'
    'Outcome: within-network FC (lag-corrected)  |  FDR-BH (35 tests)  |  '
    'Covariates: age, gender, lesion_vol, lesion_side, FD (time-varying)',
    fontsize=9, fontweight='bold')
ax_c.spines[['top','right','bottom','left']].set_visible(False)
fig_c.savefig(f'{OUT}/vph_fc_lme_heatmap.png', dpi=200, bbox_inches='tight', facecolor='white')
plt.close()
print('Saved vph_fc_lme_heatmap.png')

# ── Figure C2: Trajectory plots for top LME results ───────────────────────────
top_lme = lme_df[lme_df['p'].notna()].nsmallest(6, 'p')

if len(top_lme):
    n_panels = len(top_lme)
    ncols    = min(3, n_panels)
    nrows    = int(np.ceil(n_panels / ncols))
    fig_ct, axes_ct = plt.subplots(nrows, ncols, figsize=(5*ncols, 4.5*nrows),
                                    constrained_layout=True)
    fig_ct.suptitle('FC trajectories by oculomotor tertile (top 6 LME results, uncorrected rank)\n'
                    'Lag-corrected within-network FC', fontsize=11)
    axes_list = list(axes_ct.flat) if hasattr(axes_ct, 'flat') else [axes_ct]

    # Compute tertile cuts from acute values
    def tertile_cuts(metric):
        vals = long.loc[long.Session=='acute', f'{metric}_acute'].dropna()
        return pd.to_numeric(vals, errors='coerce').quantile([1/3, 2/3]).values

    tp_map   = {'acute': 'Acute', 'followup': '3 mo', 'followup2': '12 mo'}
    tp_order = ['Acute', '3 mo', '12 mo']

    for idx, (_, res_row) in enumerate(top_lme.iterrows()):
        ax_t   = axes_list[idx]
        metric = res_row['metric']
        net    = res_row['network']
        fc_col = f'cor_fc_{net}_{net}'
        cuts   = tertile_cuts(metric)
        color  = NET_COLORS[net]

        sub = long[['Session', f'{metric}_acute', fc_col]].copy()
        sub[f'{metric}_acute'] = pd.to_numeric(sub[f'{metric}_acute'], errors='coerce')
        sub[fc_col]            = pd.to_numeric(sub[fc_col], errors='coerce')
        sub['tertile'] = pd.cut(sub[f'{metric}_acute'],
                                bins=[-np.inf, cuts[0], cuts[1], np.inf],
                                labels=['Low','Mid','High'])
        sub['tp'] = sub['Session'].map(tp_map)

        for tertile, ls in [('Low','-'), ('Mid','--'), ('High',':')]:
            g  = sub[sub.tertile == tertile]
            mu = g.groupby('tp')[fc_col].mean().reindex(tp_order)
            se = g.groupby('tp')[fc_col].sem().reindex(tp_order)
            ax_t.plot(range(len(tp_order)), mu.values, ls=ls, color=color,
                      lw=2, marker='o', ms=6, label=tertile)
            ax_t.fill_between(range(len(tp_order)),
                               (mu-se).values, (mu+se).values,
                               color=color, alpha=0.1)

        stars_str = sig_stars(res_row.get('p_fdr', res_row['p']))
        ax_t.set_xticks(range(len(tp_order)))
        ax_t.set_xticklabels(tp_order, fontsize=8)
        ax_t.set_ylabel('Within-net FC', fontsize=8)
        ax_t.set_title(f'{METRIC_LABELS[metric]} × {NET_LABELS[net]}\n'
                       f'β={res_row["coef"]:.3f}  p_raw={res_row["p"]:.4f}{stars_str}',
                       fontsize=9, color=color, fontweight='bold')
        ax_t.legend(fontsize=7, title=METRIC_LABELS[metric], title_fontsize=7)
        ax_t.spines[['top','right']].set_visible(False)

    for idx in range(n_panels, nrows*ncols):
        axes_list[idx].set_visible(False)

    fig_ct.savefig(f'{OUT}/vph_fc_lme_trajectories.png', dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    print('Saved vph_fc_lme_trajectories.png')


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE: FC TRAJECTORY OVER RECOVERY
# 7 panels (one per within-network FC), each showing:
#   - Individual subject lines (thin, light gray)
#   - Group mean ± SEM (bold, network colour)
#   - Paired t-test significance brackets (acute→3mo, acute→12mo)
# Y-axis auto-scaled per panel. No tertile split — clean group trajectory.
# ══════════════════════════════════════════════════════════════════════════════
print('\nGenerating FC trajectory figure...')

from scipy.stats import ttest_rel

TP_SESSIONS  = ['acute', 'followup', 'followup2']
TP_LABELS    = ['Acute\n(~2 wks)', '3 months', '12 months']
TP_MAP       = dict(zip(TP_SESSIONS, TP_LABELS))
long['tp_label'] = long['Session'].map(TP_MAP)

fig_traj, axes_traj = plt.subplots(2, 4, figsize=(18, 9), constrained_layout=True)
fig_traj.suptitle(
    'Within-network FC trajectory across recovery  |  Uncorrected FC (first run)\n'
    'Bold line = group mean ± SEM  |  Gray lines = individual subjects\n'
    'Brackets: * p<0.05   ** p<0.01   *** p<0.001  (paired t-test, uncorrected)',
    fontsize=11, fontweight='bold')

axes_list_traj = list(axes_traj.flat)

def bracket(ax, x1, x2, y_top, label, color='#333333'):
    """Draw a significance bracket between x1 and x2 at height y_top."""
    h = (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.03
    ax.plot([x1, x1, x2, x2], [y_top, y_top+h, y_top+h, y_top],
            lw=1.2, color=color)
    ax.text((x1+x2)/2, y_top+h*1.2, label,
            ha='center', va='bottom', fontsize=9, color=color, fontweight='bold')

for panel_idx, net in enumerate(NETWORKS):
    ax     = axes_list_traj[panel_idx]
    fc_col = f'unc_fc_{net}_{net}'
    color  = NET_COLORS[net]

    # Build per-subject wide table
    sub = long[['ID', 'Session', fc_col]].copy()
    sub[fc_col] = pd.to_numeric(sub[fc_col], errors='coerce')

    # Individual subject trajectories (gray, thin)
    for subj_id, grp in sub.groupby('ID'):
        grp_sorted = grp.set_index('Session').reindex(TP_SESSIONS)[fc_col]
        if grp_sorted.notna().sum() >= 2:
            xpos = [i for i, s in enumerate(TP_SESSIONS) if not np.isnan(grp_sorted[s])]
            ypos = [grp_sorted[s] for s in TP_SESSIONS if not np.isnan(grp_sorted[s])]
            ax.plot(xpos, ypos, color='#aaaaaa', lw=0.5, alpha=0.35, zorder=1)

    # Group mean ± SEM
    means = [sub[sub.Session==s][fc_col].mean() for s in TP_SESSIONS]
    sems  = [sub[sub.Session==s][fc_col].sem()  for s in TP_SESSIONS]
    ns    = [sub[sub.Session==s][fc_col].notna().sum() for s in TP_SESSIONS]
    xpos  = range(3)

    ax.plot(xpos, means, color=color, lw=3, marker='o', ms=9, zorder=3)
    ax.fill_between(xpos,
                    [m-s for m,s in zip(means,sems)],
                    [m+s for m,s in zip(means,sems)],
                    color=color, alpha=0.25, zorder=2)

    # Paired t-tests (subjects with both sessions)
    def paired_t(s1, s2):
        merged = (sub[sub.Session==s1][['ID',fc_col]]
                  .merge(sub[sub.Session==s2][['ID',fc_col]], on='ID', suffixes=('_a','_b'))
                  .dropna())
        if len(merged) < 5: return 1.0
        _, p = ttest_rel(merged[f'{fc_col}_a'], merged[f'{fc_col}_b'])
        return p

    p_ac_3mo = paired_t('acute', 'followup')
    p_ac_12mo = paired_t('acute', 'followup2')

    # Auto y-top for brackets (just above max mean + SEM)
    y_top_base = max(m+s for m,s in zip(means, sems)) * 1.02
    y_rng = ax.get_ylim()

    for p_val, x1, x2, y_off in [(p_ac_3mo, 0, 1, 0.0), (p_ac_12mo, 0, 2, 0.07)]:
        if p_val < 0.05:
            star = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*'
            y_br = y_top_base + (ax.get_ylim()[1] - ax.get_ylim()[0]) * (0.08 + y_off)
            bracket(ax, x1, x2, y_br, star, color=color)

    ax.set_xticks(range(3))
    ax.set_xticklabels(TP_LABELS, fontsize=9)
    ax.set_ylabel('Within-network FC', fontsize=9)
    ax.set_title(NET_LABELS[net], fontsize=12, fontweight='bold', color=color, pad=6)

    # N labels under each timepoint
    for i, (n_val, s) in enumerate(zip(ns, TP_SESSIONS)):
        ax.text(i, ax.get_ylim()[0] - (ax.get_ylim()[1]-ax.get_ylim()[0])*0.08,
                f'n={n_val}', ha='center', fontsize=7.5, color='#666666')

    ax.spines[['top', 'right']].set_visible(False)

# Hide unused 8th panel
axes_list_traj[7].set_visible(False)

fig_traj.savefig(f'{OUT}/vph_fc_trajectory.png', dpi=200, bbox_inches='tight', facecolor='white')
plt.close()
print('Saved vph_fc_trajectory.png')


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE: Mean r(t) × FC — 7×7 heatmap across 3 timepoints
# Same format as vph_fc_unc_heatmap panel 1, repeated for acute / 3mo / 12mo
# FDR-BH across all 3 × 28 = 84 tests jointly
# ══════════════════════════════════════════════════════════════════════════════
print('\nGenerating Mean r(t) × FC trajectory heatmaps...')

SESSIONS_LONG = [('acute',      'Acute (~2 wks)'),
                 ('followup',   '3 months'),
                 ('followup2',  '12 months')]

# For each session, build a subject-level table with mean metric + FC + covariates
# Use acute lesion volume for all sessions (propagated via long dataframe)
acute_lv = (patients[patients.Session=='acute'][['ID','lesion_volume_mm3']]
            .rename(columns={'lesion_volume_mm3':'lesion_vol_acute'}))

all_records_time = []
session_dfs = {}

for sess, sess_lbl in SESSIONS_LONG:
    sess_data = patients[patients.Session == sess].copy()
    sess_data  = sess_data.merge(acute_lv, on='ID', how='left')

    for col in ['age','gender','lesion_side','lesion_vol_acute','fd','mean']:
        sess_data[col] = pd.to_numeric(sess_data[col], errors='coerce')

    covs_mean = [sess_data[c].values for c in ['age','gender','lesion_side','lesion_vol_acute','fd']
                 if c in sess_data.columns]

    for pair in FC_PAIRS:
        col = f'unc_{pair}'
        r, p, n = partial_spearman(
            sess_data['mean'].values,
            pd.to_numeric(sess_data[col], errors='coerce').values if col in sess_data.columns
            else np.full(len(sess_data), np.nan),
            covs_mean)
        parts = pair.replace('fc_','').split('_', 1)
        all_records_time.append({
            'session': sess, 'session_lbl': sess_lbl,
            'pair': pair, 'n1': parts[0], 'n2': parts[1],
            'r': r, 'p': p, 'n': n,
        })

time_df = pd.DataFrame(all_records_time)

# Joint FDR across all 84 tests
valid_t = time_df['p'].notna()
_, p_fdr_t, _, _ = multipletests(time_df.loc[valid_t,'p'], method='fdr_bh')
time_df.loc[valid_t,'p_fdr'] = p_fdr_t
time_df['stars'] = time_df['p_fdr'].apply(sig_stars)

n_sig_t = (time_df['p_fdr'] < 0.05).sum()
print(f'FDR<0.05: {n_sig_t}  |  Trends: {((time_df["p_fdr"]>=0.05)&(time_df["p_fdr"]<0.10)).sum()}  |  Total: {valid_t.sum()} tests')

# ── Build 3-panel figure ───────────────────────────────────────────────────────
fig_time, axes_time = plt.subplots(1, 3, figsize=(15, 5.5), constrained_layout=True)
fig_time.suptitle(
    'Mean r(t) × Uncorrected FC  |  Partial Spearman r  |  FDR-BH (84 tests jointly)\n'
    'Covariates: age, gender, lesion_vol, lesion_side, mean FD\n'
    'Annotated cells: FDR p<0.05 (*) or p<0.10 (†)',
    fontsize=10, fontweight='bold')

# Shared colour scale across all sessions
vmax_t = max(0.25, time_df['r'].abs().max())

for ax_t, (sess, sess_lbl) in zip(axes_time, SESSIONS_LONG):
    r7 = np.full((7,7), np.nan)
    p7 = np.full((7,7), np.nan)

    sess_rows = time_df[time_df.session == sess]
    for i, n1 in enumerate(NETWORKS):
        for j, n2 in enumerate(NETWORKS):
            pair = f'fc_{n1}_{n2}' if j >= i else f'fc_{n2}_{n1}'
            row  = sess_rows[sess_rows.pair == pair]
            if len(row):
                r7[i,j] = row.iloc[0]['r']
                r7[j,i] = row.iloc[0]['r']
                p7[i,j] = row.iloc[0].get('p_fdr', np.nan)
                p7[j,i] = row.iloc[0].get('p_fdr', np.nan)

    im = ax_t.imshow(r7, cmap='RdBu_r', vmin=-vmax_t, vmax=vmax_t, aspect='auto')
    plt.colorbar(im, ax=ax_t, fraction=0.046, label='r')

    ax_t.set_xticks(range(7)); ax_t.set_yticks(range(7))
    ax_t.set_xticklabels([NET_LABELS[n] for n in NETWORKS], fontsize=8, rotation=40, ha='right')
    ax_t.set_yticklabels([NET_LABELS[n] for n in NETWORKS], fontsize=8)
    ax_t.set_title(sess_lbl, fontsize=12, fontweight='bold')

    # Annotate only sig/trending cells
    for i in range(7):
        for j in range(7):
            if np.isnan(r7[i,j]): continue
            star = sig_stars(p7[i,j])
            if not star: continue
            color = 'white' if abs(r7[i,j]) > vmax_t*0.6 else 'black'
            ax_t.text(j, i, f'{r7[i,j]:.2f}{star}',
                      ha='center', va='center', fontsize=7.5,
                      color=color, fontweight='bold')

    # Network colour chips on y-axis
    for k, net in enumerate(NETWORKS):
        ax_t.add_patch(plt.Rectangle((-0.5, k-0.5), 0.08, 1, color=NET_COLORS[net],
                       zorder=3, transform=ax_t.get_yaxis_transform(), clip_on=False))

    ax_t.spines[['top','right','bottom','left']].set_visible(False)

fig_time.savefig(f'{OUT}/vph_fc_mean_trajectory_heatmap.png', dpi=200, bbox_inches='tight', facecolor='white')
plt.close()
print('Saved vph_fc_mean_trajectory_heatmap.png')

print('\nAll done.')

