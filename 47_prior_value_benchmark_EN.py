# -*- coding: utf-8 -*-
"""
A2: Incremental value of the senescence prior — three-arm controlled comparison.
Same GLM pipeline and resampling as the locked analysis (StratifiedKFold(5, shuffle, rs=42),
StandardScaler + LogisticRegression), 41 modelled biopsies (Norm 25 + i-IFTA 16).
  Arm C (senescence prior): the 21 ARDEGs -> 5-fold CV AUC (expect 0.896).
  Arm A (whole genome): top-k genes by |t| INSIDE each CV fold (custom transformer,
          no selection leakage) for k in {5,10,21,50,100,500,2000}.
  Arm B (random null): 100 draws of 21 genes from all expressed genes -> CV AUC distribution;
          report mean/SD/percentile of the ARDEG result within the null.
Figure: EN, 2 panels (top-k curve + random null histogram).
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
from daimon_runtime import setup_plot
setup_plot()

OUT = "results/fibrosis_ARG"
RS = 42

expr = pd.read_csv(f"{OUT}/39_GSE22459_log2expr.csv", index_col=0)
grp = pd.read_csv(f"{OUT}/39_GSE22459_groups.csv", index_col=0).iloc[:, 0]
model_idx = list(grp[grp != "IFTA"].index)
X_all = expr[model_idx].T          # 41 x G
y = (grp[model_idx] == "IFTA+I").astype(int).values
print(f"modelling set: {X_all.shape[0]} samples, {X_all.shape[1]} genes; positives={y.sum()}")

AR = pd.read_csv(f"{OUT}/39_ARDEGs.csv", index_col=0)
ARDEGS = list(AR.index)
print(f"ARDEGs: {len(ARDEGS)}")

def _rep_auc(pipe, X, y):
    """与锁定管线(脚本40)完全一致：rep 0..9, 每次 pooled OOF predict_proba -> AUC"""
    aucs = []
    for rep in range(10):
        cv_r = StratifiedKFold(5, shuffle=True, random_state=rep)
        p = cross_val_predict(pipe, X, y, cv=cv_r, method="predict_proba")[:, 1]
        aucs.append(roc_auc_score(y, p))
    return float(np.mean(aucs)), float(np.std(aucs))

def glm_cv(X, y):
    pipe = Pipeline([("sc", StandardScaler()), ("lr", LogisticRegression(max_iter=2000))])
    return _rep_auc(pipe, X, y)

class TopKSelector(BaseEstimator, TransformerMixin):
    def __init__(self, k=21):
        self.k = k
    def fit(self, X, y):
        X = np.asarray(X); y = np.asarray(y)
        t = np.abs(stats.ttest_ind(X[y == 1], X[y == 0], axis=0, equal_var=False).statistic)
        t = np.nan_to_num(t)
        self.ix_ = np.argsort(-t)[: self.k]
        return self
    def transform(self, X):
        return np.asarray(X)[:, self.ix_]

def topk_cv(X, y, k):
    pipe = Pipeline([("sel", TopKSelector(k)), ("sc", StandardScaler()),
                     ("lr", LogisticRegression(max_iter=2000))])
    return _rep_auc(pipe, X, y)

# Arm C: senescence prior (21 ARDEGs)
X_ar = X_all[ARDEGS].values
auc_ar, sd_ar = glm_cv(X_ar, y)
print(f"Arm C senescence prior (21 ARDEGs): 10x5-fold CV AUC = {auc_ar:.3f} ± {sd_ar:.3f}")

# Arm A: whole genome top-k (selection inside folds)
ks = [5, 10, 21, 50, 100, 500, 2000]
auc_topk = {}
for k in ks:
    auc_topk[k], _ = topk_cv(X_all.values, y, k)
    print(f"Arm A top-{k:>4}: 10x5-fold CV AUC = {auc_topk[k]:.3f}")

# Arm B: random 21-gene null (100 draws)
rng = np.random.RandomState(RS)
G = X_all.shape[1]
gen_arr = X_all.values
null = []
for i in range(100):
    ix = rng.choice(G, size=21, replace=False)
    null.append(glm_cv(gen_arr[:, ix], y)[0])
null = np.array(null)
pct = (null < auc_ar).mean() * 100
print(f"Arm B random-21 null: mean={null.mean():.3f}, SD={null.std():.3f}, "
      f"max={null.max():.3f}; ARDEG AUC percentile = {pct:.0f}%")

pd.DataFrame({"k": ks, "cv_auc": [auc_topk[k] for k in ks]}).to_csv(
    f"{OUT}/47_whole_genome_topk.csv", index=False, encoding="utf-8-sig")
pd.DataFrame({"random_cv_auc": null}).to_csv(
    f"{OUT}/47_random_null_21genes.csv", index=False, encoding="utf-8-sig")
with open(f"{OUT}/47_prior_benchmark_summary.csv", "w", encoding="utf-8-sig") as f:
    f.write(f"arm,setting,cv_auc\n")
    f.write(f"C,senescence_prior_21ARDEGs,{auc_ar:.4f}\n")
    for k in ks:
        f.write(f"A,whole_genome_top{k},{auc_topk[k]:.4f}\n")
    f.write(f"B,random21_mean,{null.mean():.4f}\nB,random21_sd,{null.std():.4f}\nB,random21_max,{null.max():.4f}\n")

# ---- English figure ----
fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
ax = axes[0]
ax.plot(ks, [auc_topk[k] for k in ks], "o-", color="#4a7c96", lw=1.8, label="Whole genome, top-k genes (selection inside CV folds)")
ax.axhline(auc_ar, color="#c0392b", ls="--", lw=1.8, label=f"Senescence prior (21 ARDEGs), AUC = {auc_ar:.3f}")
ax.axhline(null.mean(), color="#888888", ls=":", lw=1.2, label=f"Random 21-gene null (mean {null.mean():.2f})")
ax.axhspan(null.mean() - 2 * null.std(), null.mean() + 2 * null.std(), color="#888888", alpha=0.12, label="Random null ±2 SD")
ax.set_xscale("log")
ax.set_xlabel("Number of genes (k, log scale)")
ax.set_ylabel("5-fold cross-validation AUC")
ax.set_title("a  Discrimination vs. number of features")
ax.legend(fontsize=8.5, loc="lower right"); ax.grid(alpha=0.3)

ax = axes[1]
ax.hist(null, bins=18, color="#9db8c7", edgecolor="white", label="100 random 21-gene sets")
ax.axvline(auc_ar, color="#c0392b", lw=2, label=f"Senescence prior: {auc_ar:.3f}\n(percentile {pct:.0f}%)")
ax.set_xlabel("5-fold cross-validation AUC")
ax.set_ylabel("Count")
ax.set_title("b  Random gene-set null distribution")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(f"{OUT}/EN_figureB14_prior_value.png", dpi=150, bbox_inches="tight")
print("EN_figureB14 saved")
