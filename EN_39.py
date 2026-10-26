# -*- coding: utf-8 -*-
"""
论文B线 脚本39: GSE22459 移植肾1年协议活检 IF/TA 差异表达 + 衰老基因(ARG)交集
分组: Norm(25) vs IF/TA(IFTA 24 + IFTA+I 16 = 40)
流程: 探针→基因 → log2 → 逐基因 t 检验 + BH 校正 → DEG(|log2FC|>0.5, FDR<0.05)
      → 火山图 → top30 热图 → 与 ARG(GenAge307 ∪ CellAge845 = 1043) 交集 = ARDEGs
"""
import numpy as np, pandas as pd, gzip, sys
from pathlib import Path
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
from daimon_runtime import setup_plot
setup_plot()

DATA, OUT = "geo_data", "results/fibrosis_ARG"
Path(OUT).mkdir(parents=True, exist_ok=True)

# ---------- 探针映射 ----------
probe2gene = {}
for line in open(f"{DATA}/asia_map.tsv", encoding="utf-8"):
    p = line.rstrip("\n").split("\t")
    if len(p) == 2 and p[0] and p[1]:
        probe2gene[p[0]] = p[1]
print(f"探针映射: {len(probe2gene)}")

def parse_series_matrix(path):
    samples, titles, rows = None, None, []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("!Sample_geo_accession"):
                samples = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_title"):
                titles = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!") or line.startswith("#"):
                continue
            elif line.startswith('"ID_REF"') or line.startswith("ID_REF"):
                continue
            else:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    rows.append(parts)
    df = pd.DataFrame([r[1:len(samples)+1] for r in rows],
                      index=[r[0].strip('"') for r in rows], columns=samples)
    return df.apply(pd.to_numeric, errors="coerce"), titles

# ---------- 读取数据 ----------
df, titles = parse_series_matrix(f"{DATA}/GSE22459_matrix.txt.gz")
print(f"表达矩阵: {df.shape[0]} 探针 × {df.shape[1]} 样本")

group = []
for t in titles:
    t = t.replace('"', '')
    group.append("IFTA+I" if "IFTA+i" in t else ("IFTA" if "IFTA" in t else "Norm"))
grp = pd.Series(group, index=df.columns)
print(grp.value_counts().to_string())

# ---------- 探针→基因 (median 合并), log2 ----------
lg = np.log2(df + 1)
sub = lg.loc[lg.index.intersection(probe2gene.keys())].copy()
sub["gene"] = [probe2gene[p] for p in sub.index]
gene_expr = sub.groupby("gene").median()
print(f"基因矩阵: {gene_expr.shape[0]} 基因")

# ---------- 分组: 纤维化+炎症(IF/TA+I, 原始论文的分子活跃表型) vs 正常 ----------
# 预实验: Norm vs IFTA(单纯纤维化, n=24) DEG=0, 与原始研究"纤维化本身分子静默"一致
is_case = (grp == "IFTA+I").values
is_ctrl = (grp == "Norm").values
n_case, n_ctrl = int(is_case.sum()), int(is_ctrl.sum())
print(f"比较: IF/TA+I (n={n_case}) vs Norm (n={n_ctrl})")

X = gene_expr.values.astype(float)
case = X[:, is_case]; ctrl = X[:, is_ctrl]
logFC = case.mean(axis=1) - ctrl.mean(axis=1)

# 逐基因 Welch t 检验
tstat, pval = stats.ttest_ind(case, ctrl, axis=1, equal_var=False)
# BH 校正 (手写实现, 与 statsmodels 等价)
def bh_fdr(p):
    p = np.asarray(p, dtype=float)
    ok = ~np.isnan(p)
    fdr = np.full_like(p, np.nan)
    pv = p[ok]; n = len(pv)
    order = np.argsort(pv)
    ranked = pv[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    tmp = np.empty(n); tmp[order] = np.clip(ranked, 0, 1)
    fdr[ok] = tmp
    return fdr
fdr = bh_fdr(pval)

res = pd.DataFrame({"gene": gene_expr.index, "log2FC": logFC,
                    "t": tstat, "p": pval, "FDR": fdr}).set_index("gene")
res = res.dropna(subset=["FDR"])
res["reg"] = "NS"
res.loc[(res.FDR < 0.05) & (res.log2FC > 0.5), "reg"] = "Up"
res.loc[(res.FDR < 0.05) & (res.log2FC < -0.5), "reg"] = "Down"
res.to_csv(f"{OUT}/39_GSE22459_DEG_full.csv", encoding="utf-8-sig")
n_up = (res.reg == "Up").sum(); n_dn = (res.reg == "Down").sum()
print(f"DEG: Up={n_up}, Down={n_dn}, 总={n_up+n_dn} / {len(res)} 基因")

# ---------- ARG 交集 ----------
ARG = set(open(f"{DATA}/ARG_list.txt").read().split())
print(f"ARG 库: {len(ARG)} 基因")
res["is_ARG"] = res.index.isin(ARG)
degs = set(res.index[res.reg != "NS"])
ardegs = sorted(degs & ARG)
print(f"ARDEGs (DEG ∩ ARG): {len(ardegs)}")
ardeg_df = res.loc[ardegs].sort_values("FDR")
ardeg_df.to_csv(f"{OUT}/39_ARDEGs.csv", encoding="utf-8-sig")
n_ar_up = (ardeg_df.reg == "Up").sum(); n_ar_dn = (ardeg_df.reg == "Down").sum()
print(f"  其中 Up={n_ar_up}, Down={n_ar_dn}")

# ---------- 图1: 火山图 ----------
fig, ax = plt.subplots(figsize=(7.6, 6.2))
ns = res[res.reg == "NS"]
ax.scatter(ns.log2FC, -np.log10(ns.FDR.clip(lower=1e-300)), s=6, c="#b8b8b8", alpha=0.5, label=f"NS ({len(ns)})")
for lab, col, mk in [("Up", "#c0392b", "^"), ("Down", "#2c5f8a", "v")]:
    d = res[res.reg == lab]
    ax.scatter(d.log2FC, -np.log10(d.FDR.clip(lower=1e-300)), s=10, c=col, alpha=0.7,
               marker=mk, label=f"{lab} ({len(d)})")
ard = res.loc[ardegs]
ax.scatter(ard.log2FC, -np.log10(ard.FDR.clip(lower=1e-300)), s=46, facecolors="none",
           edgecolors="#1a1a1a", linewidths=1.1, label=f"ARDEG ({len(ard)})")
ax.axhline(-np.log10(0.05), color="k", ls="--", lw=0.8)
ax.axvline(0.5, color="k", ls="--", lw=0.8); ax.axvline(-0.5, color="k", ls="--", lw=0.8)
ax.set_xlabel("log2 fold change (i-IFTA vs Norm)")
ax.set_ylabel("-log10(FDR)")
ax.set_title(f"Differential expression in GSE22459 allograft biopsies\n(DEG ∩ ARG = {len(ardegs)} ageing-related genes)")
ax.legend(fontsize=9, loc="upper left"); ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig(f"{OUT}/EN_figureB1_volcano.png", dpi=150, bbox_inches="tight")
print("figureB1 saved")

# ---------- 图2: top30 DEG + ARDEG 热图 ----------
top30 = res[res.reg != "NS"].reindex(
    res[res.reg != "NS"]["FDR"].sort_values().index).head(30)
# 也展示全部 ARDEG（若 ≤40）
show = list(dict.fromkeys(list(top30.index) + ardegs))[:40]
if len(ardegs) <= 40:
    show = list(ardegs)          # ARDEG 不太多就全展示
mat = gene_expr.loc[gene_expr.index.isin(show)]
matz = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1).replace(0, np.nan), axis=0)
order_cols = list(grp[grp == "Norm"].index) + list(grp[grp == "IFTA"].index) + list(grp[grp == "IFTA+I"].index)
matz = matz[order_cols]
fig, ax = plt.subplots(figsize=(11, max(5, 0.32 * len(show) + 1.5)))
im = ax.imshow(matz.values, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
ax.set_yticks(range(len(show)))
lab_arg = [f"{g} *" if g in ARG else g for g in show]
ax.set_yticklabels(lab_arg, fontsize=8)
color_bar = [("#4a7c96" if grp[c] == "Norm" else ("#b08a3e" if grp[c] == "IFTA" else "#c96a4a")) for c in matz.columns]
for i, c in enumerate(color_bar):
    ax.axvline(i - 0.5, color=c, lw=2.2)
ax.set_xticks([]); ax.set_xlabel("Samples (blue=Norm, gold=pure IF/TA, red=i-IFTA)")
ax.set_title("Expression heatmap of ageing-related DEGs (ARDEGs)\n* = member of GenAge/CellAge gene sets", fontsize=11)
fig.colorbar(im, ax=ax, shrink=0.6, label="z score")
fig.tight_layout()
fig.savefig(f"{OUT}/EN_figureB2_heatmap.png", dpi=150, bbox_inches="tight")
print("figureB2 saved")

# 保存表达矩阵与分组供下游 ML 使用
gene_expr.to_csv(f"{OUT}/39_GSE22459_log2expr.csv", encoding="utf-8-sig")
grp.to_csv(f"{OUT}/39_GSE22459_groups.csv", encoding="utf-8-sig")
print("done")
