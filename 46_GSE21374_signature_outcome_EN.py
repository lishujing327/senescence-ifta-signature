# -*- coding: utf-8 -*-
"""
A1: GSE21374 outcome linkage — 5-gene senescence signature vs graft failure.
Cohort: late (>1 yr) first biopsy per patient (late_main==1), GPL570.
Signature: locked coefficients from 40_signature_coefficients.csv, cohort-wise z-weighted.
Stats: Cox PH (Breslow, hand-written as in script 10), KM by median split, log-rank.
English figure output for submission.
"""
import numpy as np, pandas as pd, gzip, sys, csv
from pathlib import Path
from scipy import stats, optimize
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
from daimon_runtime import setup_plot
setup_plot()

DATA, OUT = "geo_data", "results/fibrosis_ARG"

# ---- locked signature ----
coef = {}
with open(f"{OUT}/40_signature_coefficients.csv", encoding="utf-8-sig") as f:
    for row in csv.DictReader(f):
        coef[row["gene"]] = float(row["beta"])
GENES = list(coef)
print("signature genes:", GENES)

probe2gene = {}
for line in open(f"{DATA}/asia_map.tsv", encoding="utf-8"):
    p = line.rstrip("\n").split("\t")
    if len(p) >= 2 and p[0] and p[1]:
        probe2gene[p[0].strip('"')] = p[1].strip('"')

def parse_series(path):
    samples, chars, rows = None, [], []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("!Sample_geo_accession"):
                samples = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                chars.append([x.strip('"') for x in line.rstrip("\n").split("\t")[1:]])
            elif line.startswith("!") or line.startswith("#") or line.startswith('"ID_REF"') or line.startswith("ID_REF"):
                continue
            else:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    rows.append(parts)
    df = pd.DataFrame([r[1:len(samples)+1] for r in rows],
                      index=[r[0].strip('"') for r in rows], columns=samples)
    df = df.apply(pd.to_numeric, errors="coerce")
    meta = pd.DataFrame(index=samples)
    for ch in chars:
        for i, v in enumerate(ch):
            if ": " in v:
                k, val = v.split(": ", 1)
                meta.loc[samples[i], k] = val
    return df, meta

df, meta = parse_series(f"{DATA}/GSE21374_series_matrix.txt.gz")
print(f"matrix {df.shape}")

clin = pd.DataFrame(index=df.columns)
clin["failed"] = pd.to_numeric(meta["failed=1/non failed=0"], errors="coerce")
clin["time"] = pd.to_numeric(meta["time from biopsy to failure/censoring (days)"], errors="coerce")
clin["late_main"] = pd.to_numeric(
    meta["sample included in the main analysis- the first biopsy per patient from the 105 late (> 1 year) patients"],
    errors="coerce")

# ---- signature score ----
expr = np.log2(df + 1)
sub = expr.loc[expr.index.isin(probe2gene.keys())].copy()
sub["gene"] = [probe2gene[p] for p in sub.index]
genes = sub.groupby("gene").mean()
present = [g for g in GENES if g in genes.index]
print(f"genes mapped: {len(present)}/{len(GENES)} -> {present}")
z = genes.sub(genes.mean(axis=1), axis=0).div(genes.std(axis=1).replace(0, np.nan), axis=0)
sig = sum(coef[g] * z.loc[g] for g in present)
clin["sig"] = sig

main = clin[(clin["late_main"] == 1) & clin["time"].notna() & clin["failed"].notna()].copy()
print(f"main cohort n={len(main)}, events={int(main['failed'].sum())}, "
      f"median follow-up {main['time'].median():.0f} d")

# ---- hand-written survival ----
def km_curve(time, event):
    t = np.asarray(time, float); e = np.asarray(event, float)
    ts = np.sort(np.unique(t)); surv = []; at_risk = len(t); S = 1.0
    for x in ts:
        d = ((t == x) & (e == 1)).sum(); c = ((t == x) & (e == 0)).sum()
        if at_risk > 0:
            S *= (1 - d / at_risk)
        surv.append(S); at_risk -= (d + c)
    return ts, np.array(surv)

def logrank(time, event, group):
    t = np.asarray(time, float); e = np.asarray(event, float); g = np.asarray(group)
    ts = np.sort(np.unique(t[e == 1])); O1 = E1 = V = 0.0
    for x in ts:
        at = t >= x; n = at.sum(); n1 = (at & (g == 1)).sum(); n0 = n - n1
        d = ((t == x) & (e == 1)).sum()
        if n <= 1: continue
        E1 += d * n1 / n
        V += d * (n1 / n) * (n0 / n) * (n - d) / (n - 1)
        O1 += ((t == x) & (e == 1) & (g == 1)).sum()
    z2 = (O1 - E1) ** 2 / V if V > 0 else np.nan
    from scipy.stats import chi2
    return z2, 1 - chi2.cdf(z2, 1)

def cox(time, event, x):
    t = np.asarray(time, float); e = np.asarray(event, float); X = np.asarray(x, float)
    X = (X - X.mean()) / X.std()
    order = np.argsort(t); t, e, X = t[order], e[order], X[order]
    def nll(beta):
        lp = X * beta[0]
        ex = np.exp(lp); r = np.cumsum(ex[::-1])[::-1]
        ev = e == 1
        ll = (lp[ev]).sum() - (np.log(r[ev])).sum()
        return -ll
    res = optimize.minimize(nll, [0.0], method="BFGS")
    beta = res.x[0]
    eps = 1e-5
    n2 = nll([beta + eps]); n0v = nll([beta - eps])
    se = np.sqrt(eps ** 2 / max(n2 - 2 * nll([beta]) + n0v, 1e-12))
    from scipy.stats import norm
    z = beta / se; p = 2 * (1 - norm.cdf(abs(z)))
    lo = (beta - 1.96 * se); hi = (beta + 1.96 * se)
    return np.exp(beta), np.exp(lo), np.exp(hi), p

# per 1-SD HR
sd = main["sig"].std()
hr, lo, hi, p = cox(main["time"], main["failed"], main["sig"])
print(f"Cox per 1 SD: HR={hr:.2f} (95%CI {lo:.2f}-{hi:.2f}), p={p:.4f}")

# KM by median
med = main["sig"].median()
main["grp"] = (main["sig"] >= med).astype(int)
lr_z2, lr_p = logrank(main["time"], main["failed"], main["grp"])
print(f"log-rank (median split): chi2={lr_z2:.2f}, p={lr_p:.4f}")

# save
main.to_csv(f"{OUT}/46_GSE21374_signature_outcome.csv", encoding="utf-8-sig")
with open(f"{OUT}/46_GSE21374_HR_summary.csv", "w", encoding="utf-8-sig") as f:
    f.write("cohort,n,events,HR_per_SD,HR_lo,HR_hi,p,logrank_chi2,logrank_p\n")
    f.write(f"GSE21374 late-first-biopsy,{len(main)},{int(main['failed'].sum())},"
            f"{hr:.3f},{lo:.3f},{hi:.3f},{p:.5f},{lr_z2:.3f},{lr_p:.6f}\n")

# ---- English figure ----
fig, axes = plt.subplots(1, 2, figsize=(12.5, 5))
ax = axes[0]
for g, lab, c in [(0, f"Low signature (n={(main['grp']==0).sum()})", "#4a7c96"),
                  (1, f"High signature (n={(main['grp']==1).sum()})", "#c96a4a")]:
    m = main[main["grp"] == g]
    ts, sv = km_curve(m["time"], m["failed"])
    years = ts / 365.25
    ax.step(years, sv, where="post", color=c, lw=1.8, label=lab)
    ax.plot(years[-1], sv[-1], "o", color=c, ms=4)
ax.set_xlabel("Years from biopsy")
ax.set_ylabel("Graft survival probability")
lr_p_txt = "P < 0.001" if lr_p < 0.001 else f"P = {lr_p:.3f}"
ax.set_title(f"a  Kaplan-Meier by median signature score\n(log-rank {lr_p_txt})")
ax.legend(fontsize=9, loc="lower left"); ax.grid(alpha=0.3)
ax.set_ylim(0, 1.02)

ax = axes[1]
ax.hist(main.loc[main["failed"] == 0, "sig"], bins=15, alpha=0.6, color="#4a7c96", label="No graft failure")
ax.hist(main.loc[main["failed"] == 1, "sig"], bins=8, alpha=0.7, color="#c0392b", label="Graft failure")
ax.axvline(med, color="k", ls="--", lw=1)
ax.set_xlabel("Senescence signature score (cohort-wise z-weighted)")
ax.set_ylabel("Number of biopsies")
ax.set_title(f"b  Signature score distribution\n(Cox HR per 1 SD = {hr:.2f}, P = {p:.3f})")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(f"{OUT}/EN_figureB13_outcome_GSE21374.png", dpi=150, bbox_inches="tight")
print("EN_figureB13 saved")
