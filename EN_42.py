# -*- coding: utf-8 -*-
"""
42_外部验证.py
用锁定的5基因衰老签名（40_signature_coefficients.csv）在两个独立外部队列验证：
  GSE53605 (GPL571, n=55): IF/TA 10 / Normal 18 / ACR 13 / CNIT 14 —— 主验证 + 特异性
  GSE7392  (GPL570, n=30): 亚临床间质纤维化(GIF/TA) 16 / Normal 14 —— 同平台二次验证
签名得分 = sum(beta * zscore(基因表达))，队列内逐基因 z 化。
"""
import gzip, csv, sys, io
from pathlib import Path
import numpy as np
from scipy import stats
from sklearn.metrics import roc_curve, auc

sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
from daimon_runtime import setup_plot
import matplotlib
matplotlib.use("Agg")

OUT = Path("results/fibrosis_ARG"); OUT.mkdir(parents=True, exist_ok=True)

# ---------- 签名 ----------
coef = {}
with open(OUT/"40_signature_coefficients.csv", encoding="utf-8-sig") as f:
    for row in csv.DictReader(f):
        coef[row["gene"]] = float(row["beta"])
GENES = list(coef)
print("签名基因:", GENES)

# ---------- 探针映射（GPL570 注释，GPL571 探针为其子集） ----------
probe2sym = {}
with open("geo_data/asia_map.tsv", encoding="utf-8") as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2 and p[0] and p[1]:
            probe2sym[p[0].strip('"')] = p[1].strip('"')
print("映射表探针数:", len(probe2sym))

def read_series(path):
    """返回 (titles, chars(list of lists), {gene: np.array})"""
    titles, chars, table, in_tab = [], [], [], False
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                titles = [t.strip('"') for t in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                chars.append([c.strip('"') for c in line.rstrip("\n").split("\t")[1:]])
            elif line.startswith("!series_matrix_table_begin"):
                in_tab = True; continue
            elif line.startswith("!series_matrix_table_end"):
                break
            elif in_tab:
                p = line.rstrip("\n").split("\t")
                if p[0].startswith('"ID') or not p[0]:
                    continue
                table.append(p)
    n = len(titles)
    gene_expr = {}
    for p in table:
        sym = probe2sym.get(p[0].strip('"'))
        if not sym: continue
        try:
            v = np.array([float(x) if x not in ("", "null") else np.nan for x in p[1:n+1]])
        except ValueError:
            continue
        if np.all(np.isnan(v)): continue
        if sym in gene_expr:  # 多探针取均值
            gene_expr[sym] = np.nanmean(np.vstack([gene_expr[sym], v]), axis=0)
        else:
            gene_expr[sym] = v
    return titles, chars, gene_expr

def signature_score(gene_expr, do_log2):
    """逐基因 z 化后加权和；返回 (score, present_genes)"""
    X = []
    present = []
    for g in GENES:
        if g not in gene_expr:
            print(f"  !! 缺基因 {g}"); continue
        v = np.asarray(gene_expr[g], dtype=float)
        if do_log2:
            v = np.log2(np.clip(v, 0, None) + 1.0)
        m = np.nanmean(v); s = np.nanstd(v)
        if s == 0 or np.isnan(s): s = 1.0
        X.append((v - m) / s)
        present.append(g)
    return np.sum([coef[g] * x for g, x in zip(present, X)], axis=0), present

def boot_auc(y, s, n_boot=2000, seed=42):
    rng = np.random.RandomState(seed)
    aucs = []
    y = np.asarray(y); s = np.asarray(s)
    pos = s[y == 1]; neg = s[y == 0]
    for _ in range(n_boot):
        pi = rng.randint(0, len(pos), len(pos)); ni = rng.randint(0, len(neg), len(neg))
        a = auc_from_samples(pos[pi], neg[ni])
        if not np.isnan(a): aucs.append(a)
    a0 = auc_from_samples(pos, neg)
    lo, hi = np.percentile(aucs, [2.5, 97.5])
    return a0, lo, hi

def auc_from_samples(pos, neg):
    r = stats.rankdata(np.concatenate([pos, neg]))
    r1 = r[:len(pos)].sum()
    return (r1 - len(pos)*(len(pos)+1)/2) / (len(pos)*len(neg))

results = []

# ---------- GSE53605 ----------
titles, chars, ge = read_series("geo_data/GSE53605_matrix.txt.gz")
diag = chars[1] if len(chars) > 1 else chars[0]   # histological diagnosis 行
print("\nGSE53605 表型行示例:", diag[:3])
def gse53605_group(d):
    d = d.lower()
    if "ifta" in d or "if/ta" in d: return "IF/TA"
    if "normal" in d: return "Normal"
    if "acr" in d or "acute rejection" in d: return "ACR"
    if "cnit" in d: return "CNIT"
    return "?"
groups53605 = [gse53605_group(d) for d in diag]
from collections import Counter
print("GSE53605 分组:", Counter(groups53605))
score5, present = signature_score(ge, do_log2=False)
print("  命中基因:", present)
groups53605 = np.array(groups53605)
np.save(OUT/"_tmp_score5.npy", score5)
with open(OUT/"42_EN_GSE53605_scores.csv", "w", newline="", encoding="utf-8-sig") as f:
    w = csv.writer(f); w.writerow(["title", "group", "signature_score"])
    for t, g, s in zip(titles, groups53605, score5):
        w.writerow([t, g, round(float(s), 4)])

for pos_label, neg_label in [("IF/TA", "Normal"), ("ACR", "Normal"), ("CNIT", "Normal")]:
    mask = np.isin(groups53605, [pos_label, neg_label])
    y = (groups53605[mask] == pos_label).astype(int)
    s = score5[mask]
    a, lo, hi = boot_auc(y, s)
    p = stats.mannwhitneyu(s[y == 1], s[y == 0]).pvalue
    results.append(dict(cohort="GSE53605", comparison=f"{pos_label} vs Normal",
                        n_pos=int(y.sum()), n_neg=int((1-y).sum()), auc=round(a,3),
                        ci=f"{lo:.3f}-{hi:.3f}", p=p))
    print(f"GSE53605 {pos_label} vs Normal: AUC={a:.3f} ({lo:.3f}-{hi:.3f}), p={p:.2e}, n={y.sum()}/{(1-y).sum()}")

# ---------- GSE7392 ----------
titles2, chars2, ge2 = read_series("geo_data/GSE7392_matrix.txt.gz")
groups7392 = np.array(["GIF/TA" if t.startswith("Subclinical") else "Normal" for t in titles2])
print("\nGSE7392 分组:", Counter(groups7392))
score7, present7 = signature_score(ge2, do_log2=True)
print("  命中基因:", present7)
with open(OUT/"42_EN_GSE7392_scores.csv", "w", newline="", encoding="utf-8-sig") as f:
    w = csv.writer(f); w.writerow(["title", "group", "signature_score"])
    for t, g, s in zip(titles2, groups7392, score7):
        w.writerow([t, g, round(float(s), 4)])
mask = np.ones(len(groups7392), bool)
y = (groups7392 == "GIF/TA").astype(int)
a, lo, hi = boot_auc(y, score7)
p = stats.mannwhitneyu(score7[y == 1], score7[y == 0]).pvalue
results.append(dict(cohort="GSE7392", comparison="GIF/TA vs Normal",
                    n_pos=int(y.sum()), n_neg=int((1-y).sum()), auc=round(a,3),
                    ci=f"{lo:.3f}-{hi:.3f}", p=p))
print(f"GSE7392 GIF/TA vs Normal: AUC={a:.3f} ({lo:.3f}-{hi:.3f}), p={p:.2e}, n={y.sum()}/{(1-y).sum()}")

# 患者级（GSE7392 为配对设计：每例 T=0/T=12m 各一）
pat = np.array([t.split("(")[-1].split(" ")[0].rstrip("#0123456789") + "#" + t.split("#")[-1].split(" ")[0] for t in titles2])
# 更稳妥：从标题 "GIF#1 T=0" 提取患者号
import re
pat_ids = []
for t in titles2:
    m = re.search(r"(GIF|Norm)#(\d+)", t)
    pat_ids.append(f"{m.group(1)}{m.group(2)}" if m else t)
pat_ids = np.array(pat_ids)
pat_score = {pid: np.mean(score7[pat_ids == pid]) for pid in np.unique(pat_ids)}
pat_group = {pid: ("GIF/TA" if pid.startswith("GIF") else "Normal") for pid in pat_score}
ps = np.array([pat_score[k] for k in pat_score]); pg = np.array([1 if pat_group[k]=="GIF/TA" else 0 for k in pat_score])
a_p = auc_from_samples(ps[pg==1], ps[pg==0])
p_p = stats.mannwhitneyu(ps[pg==1], ps[pg==0]).pvalue
results.append(dict(cohort="GSE7392(患者级)", comparison="GIF/TA vs Normal",
                    n_pos=int(pg.sum()), n_neg=int((1-pg).sum()), auc=round(a_p,3),
                    ci="-", p=p_p))
print(f"GSE7392 患者级: AUC={a_p:.3f}, p={p_p:.2e}")

with open(OUT/"42_EN_validation_AUC.csv", "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=["cohort", "comparison", "n_pos", "n_neg", "auc", "ci", "p"])
    w.writeheader()
    for r in results:
        r2 = dict(r); r2["p"] = f"{r['p']:.2e}"; w.writerow(r2)

# ---------- 图 figureB12 ----------
setup_plot()
import matplotlib.pyplot as plt
fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))

# (a) GSE53605 四组箱线
ax = axes[0]
order = ["Normal", "IF/TA", "ACR", "CNIT"]
colors = ["#9ecae1", "#e6550d", "#bdbdbd", "#bdbdbd"]
data = [score5[groups53605 == g] for g in order]
bp = ax.boxplot(data, labels=[f"{g}\n(n={len(d)})" for g, d in zip(order, data)],
                patch_artist=True, widths=0.55, medianprops=dict(color="black", lw=1.4))
for patch, c in zip(bp["boxes"], colors):
    patch.set_facecolor(c); patch.set_alpha(0.75)
ax.axhline(0, ls="--", lw=0.8, c="gray")
ax.set_ylabel("Senescence signature score (cohort-wise z-weighted)")
ax.set_title("a  GSE53605 (GPL571, n=55)", fontsize=11)
ax.set_xlabel("Histological diagnosis")

# (b) GSE53605 ROC × 3
ax = axes[1]
for pos_label, c, lb in [("IF/TA", "#e6550d", None), ("ACR", "#756bb1", None), ("CNIT", "#31a354", None)]:
    mask = np.isin(groups53605, [pos_label, "Normal"])
    yy = (groups53605[mask] == pos_label).astype(int); ss = score5[mask]
    fpr, tpr, _ = roc_curve(yy, ss)
    a = auc(fpr, tpr)
    ax.plot(fpr, tpr, c=c, lw=1.8, label=f"{pos_label} vs N (AUC {a:.2f})")
ax.plot([0, 1], [0, 1], ls="--", c="gray", lw=0.8)
ax.set_xlabel("1 − specificity"); ax.set_ylabel("Sensitivity")
ax.set_title("b  GSE53605 signature ROC", fontsize=11)
ax.legend(fontsize=8.5, loc="lower right")

# (c) GSE7392 箱线 + 患者级均值
ax = axes[2]
order2 = ["Normal", "GIF/TA"]
data2 = [score7[groups7392 == g] for g in order2]
bp = ax.boxplot(data2, labels=[f"{g}\n(n={len(d)})" for g, d in zip(order2, data2)],
                patch_artist=True, widths=0.45, medianprops=dict(color="black", lw=1.4))
for patch, c in zip(bp["boxes"], ["#9ecae1", "#e6550d"]):
    patch.set_facecolor(c); patch.set_alpha(0.75)
# 配对患者连线（T=0 → T=12m）
for pid in np.unique(pat_ids):
    m = pat_ids == pid
    if m.sum() == 2:
        ts = np.array([t for t in np.array(titles2)[m]])
        order_t = np.argsort(ts)
        ax.plot([1, 2][::1], score7[m][order_t], c="gray", lw=0.7, alpha=0.55, zorder=1)
ax.axhline(0, ls="--", lw=0.8, c="gray")
ax.set_ylabel("Senescence signature score")
ax.set_title("c  GSE7392 (GPL570, paired biopsies, n=30)", fontsize=11)
ax.set_xlabel("Diagnosis (grey lines = same patient, T=0→12 mo)")

fig.tight_layout()
fig.savefig(OUT/"EN_figureB12_external_validation.png", dpi=200, bbox_inches="tight")
print("\nfigureB12 已保存")
print("\n=== AUC 汇总 ===")
for r in results:
    print(r)
