# -*- coding: utf-8 -*-
"""
A6: GSE147451 (192例死亡供肾植入前活检, GPL571) —— 供体年龄相关性
问题: 受体侧5基因签名/21-ARDEG衰老活性评分,在供体侧是否随供体年龄增长?
口径: 锁定签名系数不重构; 队列内z加权; Spearman相关 + 偏相关(校正batch effect)
      + >=70 / >=75 年龄组比较 + 与KDPI相关
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
import gzip
import numpy as np, pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from daimon_runtime import setup_plot
setup_plot()

OUT = "results/纤维化ARG"

# ---------- 1. 解析series matrix ----------
f = gzip.open("geo_data/GSE147451_matrix.txt.gz", "rt", encoding="utf-8", errors="replace")
gsms, ch_rows, table = None, [], None
in_table = False
for line in f:
    line = line.rstrip("\n")
    if line.startswith("!Sample_geo_accession"):
        gsms = [x.strip(' "') for x in line.split("\t")[1:]]
    elif line.startswith("!Sample_characteristics_ch1"):
        ch_rows.append([x.strip(' "') for x in line.split("\t")[1:]])
    elif line.startswith("!series_matrix_table_begin"):
        in_table = True
        header = next(f).rstrip("\n").split("\t")
        rows = []
        for tl in f:
            tl = tl.rstrip("\n")
            if tl.startswith("!series_matrix_table_end"):
                break
            rows.append(tl.split("\t"))
f.close()
expr = pd.DataFrame([[float(v) for v in r[1:]] for r in rows],
                    index=[r[0].strip('"') for r in rows], columns=[h.strip('"') for h in header[1:]])
print("matrix:", expr.shape)

# 每样本特征(按列汇总6行characteristics)
meta = {}
for j, g in enumerate(gsms):
    d = {}
    for row_ in ch_rows:
        if j < len(row_) and row_[j]:
            k, _, v = row_[j].partition(": ")
            d[k.strip()] = v.strip()
    meta[g] = d
meta = pd.DataFrame(meta).T
print("meta fields:", list(meta.columns), "| n =", len(meta))

# ---------- 2. 探针映射 + 基因表达 ----------
probe2gene = {}
for line in open("geo_data/asia_map.tsv", encoding="utf-8"):
    p = line.rstrip("\n").split("\t")
    if len(p) >= 2 and p[0]:
        probe2gene[p[0]] = p[1]
sub = expr.loc[[p for p in expr.index if p in probe2gene]]
gene_expr = sub.assign(gene=[probe2gene[p] for p in sub.index]).groupby("gene").mean()
print("genes after mapping:", gene_expr.shape[0], "| TLR10 present:", "TLR10" in gene_expr.index)

if gene_expr.values.max() > 50:
    gene_expr = np.log2(gene_expr + 1)
    print("log2 transformed")

# ---------- 3. 签名评分(锁定系数, 队列内z加权) ----------
coef = pd.read_csv(f"{OUT}/40_EN_signature_coefficients.csv", index_col=0)["beta"]
genes_5 = [g for g in coef.index if g in gene_expr.index]
print(f"signature genes available: {genes_5}")
z = gene_expr.sub(gene_expr.mean(axis=1), axis=0).div(gene_expr.std(axis=1), axis=0)
score5 = (z.loc[genes_5].mul(coef.loc[genes_5], axis=0)).sum()

# 21-ARDEG 衰老活性评分
ardegs = pd.read_csv(f"{OUT}/39_ARDEGs.csv", index_col=0).index.tolist()
ar = [g for g in ardegs if g in gene_expr.index]
score21 = z.loc[ar].mean()
print(f"21-ARDEG score uses {len(ar)} genes")

df = pd.DataFrame({"donor_age": pd.to_numeric(meta["donor age"], errors="coerce"),
                   "kdpi": pd.to_numeric(meta["kdpi"], errors="coerce"),
                   "batch": pd.to_numeric(meta["batch effect"], errors="coerce"),
                   "score5": score5, "score21": score21}).dropna(subset=["donor_age"])
df["age70"] = (df.donor_age >= 70).map({True: ">=70y", False: "<70y"})
df["age75"] = (df.donor_age >= 75).map({True: ">=75y", False: "<75y"})
print("analysis n =", len(df), "| age range:", int(df.donor_age.min()), "-", int(df.donor_age.max()))

# ---------- 4. 统计 ----------
def spear(a, b):
    r, p = stats.spearmanr(a, b)
    return r, p

def partial_spearman(x, y, cov):
    rx = stats.linregress(cov, x).resid if False else None
    # 残差法偏相关
    rx = x - np.poly1d(np.polyfit(cov, x, 1))(cov)
    ry = y - np.poly1d(np.polyfit(cov, y, 1))(cov)
    return stats.spearmanr(rx, ry)

rows = []
for name, s in [("5-gene signature", df.score5), ("21-ARDEG score", df.score21)]:
    r, p = spear(df.donor_age, s)
    sub_b = df.batch.notna()
    if sub_b.sum() > 30:
        rb, pb = partial_spearman(df.donor_age[sub_b].values, s[sub_b].values, df.batch[sub_b].values)
    else:
        rb, pb = (np.nan, np.nan)
    rk, pk = spear(df.kdpi.dropna(), s[df.kdpi.notna()]) if df.kdpi.notna().sum() > 10 else (np.nan, np.nan)
    rows.append([name, len(df), round(r, 3), p, round(rb, 3) if sub_b.sum() > 30 else np.nan, pb, round(rk, 3), pk])

for label, col in [(">=70y vs <70y", "age70"), (">=75y vs <75y", "age75")]:
    for name, s in [("5-gene signature", df.score5), ("21-ARDEG score", df.score21)]:
        a = s[df[col] == df[col].unique()[0] if False else df[col].str.startswith(">=")]
        b = s[df[col].str.startswith("<")]
        u, p = stats.mannwhitneyu(a, b)
        rows.append([f"{name} | {label}", f"{len(a)} vs {len(b)}", round(a.median(), 3), p, np.nan, np.nan, np.nan, np.nan])

res = pd.DataFrame(rows, columns=["analysis", "n", "spearman_r_or_median_high", "p", "partial_r_batchadj", "partial_p", "kdpi_spearman_r", "kdpi_p"])
res.to_csv(f"{OUT}/43_GSE147451_donor_age.csv", index=False, encoding="utf-8-sig")
print(res.to_string(index=False))
df.to_csv(f"{OUT}/43_GSE147451_scores.csv", encoding="utf-8-sig")

# ---------- 5. 图 ----------
fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
for ax, (name, s) in zip(axes, [("4-gene signature score (TLR10 probe absent on this platform)", df.score5), ("21-ARDEG senescence-activity score", df.score21)]):
    young = df.donor_age < 70; old = df.donor_age >= 70
    ax.scatter(df.donor_age[young], s[young], c="#4a7c96", s=28, alpha=0.75, label="<70y")
    ax.scatter(df.donor_age[old], s[old], c="#b04a4a", s=34, alpha=0.9, label=">=70y")
    k = np.polyfit(df.donor_age, s, 1)
    xs = np.linspace(df.donor_age.min(), df.donor_age.max(), 50)
    ax.plot(xs, np.poly1d(k)(xs), color="#555555", lw=1.4, ls="--")
    r, p = stats.spearmanr(df.donor_age, s)
    ax.set_title(f"{name}\nSpearman r = {r:.3f}, P = {p:.3g}", fontsize=10.5)
    ax.set_xlabel("Donor age (years)"); ax.grid(alpha=0.3); ax.legend(fontsize=8.5, framealpha=0.7)
axes[0].set_ylabel("Locked signature score (cohort z-weighted)")
fig.suptitle("GSE147451: donor-side pre-implantation biopsies (n = 192, deceased donors)", fontsize=11.5, y=1.02)
fig.tight_layout()
fig.savefig(f"{OUT}/EN_figureB17_donor_age_GSE147451.png", dpi=300, bbox_inches="tight")
print("saved: EN_figureB17_donor_age_GSE147451.png, 43_GSE147451_donor_age.csv, 43_GSE147451_scores.csv")
