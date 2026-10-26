# -*- coding: utf-8 -*-
"""
论文B线 脚本40: 21个ARDEGs诊断IF/TA+I的机器学习建模
样本: GSE22459 中 Norm(25) vs IFTA+I(16), IFTA(24)仅用于评分谱展示不入模
模型: 逻辑回归(GLM) / SVM-RBF / 随机森林 / 梯度提升  (sklearn实现; xgboost未安装, GBM替代并注明)
评估: 重复5折CV AUC, 特征重要性共识排名 → top5枢纽基因
输出: 模型AUC对比图 + 单基因AUC + 组合签名ROC + 列线图 + 校准曲线 + DCA + 三组评分谱
"""
import numpy as np, pandas as pd, sys
from pathlib import Path
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_curve, auc, roc_auc_score
sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
from daimon_runtime import setup_plot
setup_plot()

OUT = "results/fibrosis_ARG"
RNG = np.random.RandomState(42)

expr = pd.read_csv(f"{OUT}/39_GSE22459_log2expr.csv", index_col=0)
grp = pd.read_csv(f"{OUT}/39_GSE22459_groups.csv", index_col=0).iloc[:, 0]
genes = pd.read_csv(f"{OUT}/39_ARDEGs.csv", index_col=0).index.tolist()
print(f"建模基因 {len(genes)}; 样本组: {grp.value_counts().to_dict()}")

train_mask = grp.isin(["Norm", "IFTA+I"])
train_cols = grp[train_mask].index                       # 41个入模样本
X = expr.loc[genes, train_cols].T.values                # (41, 21)
y = (grp[train_mask] == "IFTA+I").astype(int).values
names = np.array(genes)
print(f"训练: n={len(y)}, 阳性={y.sum()}")

# ---------- 四模型 重复5折CV ----------
models = {
    "GLM (logistic)": make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000)),
    "SVM":           make_pipeline(StandardScaler(), SVC(probability=True, random_state=42)),
    "Random forest":  RandomForestClassifier(n_estimators=500, random_state=42),
    "GBM":            GradientBoostingClassifier(random_state=42),
}
cv = StratifiedKFold(5, shuffle=True, random_state=42)
rows = []
for mname, mdl in models.items():
    aucs = []
    for rep in range(10):
        cv_r = StratifiedKFold(5, shuffle=True, random_state=rep)
        p = cross_val_predict(mdl, X, y, cv=cv_r, method="predict_proba")[:, 1]
        aucs.append(roc_auc_score(y, p))
    rows.append((mname, np.mean(aucs), np.std(aucs)))
    print(f"{mname}: CV AUC = {np.mean(aucs):.3f} ± {np.std(aucs):.3f}")
res = pd.DataFrame(rows, columns=["model", "auc_mean", "auc_sd"]).sort_values("auc_mean", ascending=False)
res.to_csv(f"{OUT}/40_EN_model_AUC.csv", index=False, encoding="utf-8-sig")

# ---------- 特征重要性共识排名 ----------
imp = pd.DataFrame(index=names)
rf = RandomForestClassifier(n_estimators=1000, random_state=42).fit(X, y)
imp["RF"] = rf.feature_importances_
lr = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000)).fit(X, y)
imp["GLM"] = np.abs(lr.named_steps["logisticregression"].coef_[0])
gb = GradientBoostingClassifier(random_state=42).fit(X, y)
imp["GBM"] = gb.feature_importances_
for c in imp.columns:
    imp[c + "_rank"] = imp[c].rank(ascending=False)
imp["mean_rank"] = imp[[c + "_rank" for c in ["RF", "GLM", "GBM"]]].mean(axis=1)
imp = imp.sort_values("mean_rank")
imp.to_csv(f"{OUT}/40_EN_feature_importance.csv", encoding="utf-8-sig")
HUB = imp.head(5).index.tolist()
print("top5 枢纽基因:", HUB)

# ---------- 组合签名 (GLM on hub genes) ----------
Xh = expr.loc[HUB, train_cols].T.values
sig = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000)).fit(Xh, y)
p_all = sig.predict_proba(Xh)[:, 1]
cv_p = cross_val_predict(sig, Xh, y, cv=cv, method="predict_proba")[:, 1]
cv_auc = roc_auc_score(y, cv_p)
fpr_t, tpr_t, _ = roc_curve(y, p_all)
fpr_c, tpr_c, _ = roc_curve(y, cv_p)
print(f"5基因签名: 训练AUC={auc(fpr_t, tpr_t):.3f}, 5折CV AUC={cv_auc:.3f}")

# 单基因 AUC
single = []
for g in HUB:
    xg = expr.loc[g, train_cols].values.reshape(-1, 1)
    pg = cross_val_predict(make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000)),
                           xg, y, cv=cv, method="predict_proba")[:, 1]
    single.append((g, roc_auc_score(y, pg)))
single = pd.DataFrame(single, columns=["gene", "cv_auc"]).sort_values("cv_auc", ascending=False)
single.to_csv(f"{OUT}/40_EN_single_gene_AUC.csv", index=False, encoding="utf-8-sig")
print(single.to_string(index=False))

# ---------- 图1: 模型AUC对比 + 特征重要性 ----------
fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
ax = axes[0]
cols = ["#4a7c96", "#b08a3e", "#7c8f4a", "#8A5A44"]
ax.bar(res.model, res.auc_mean, yerr=res.auc_sd, color=cols, alpha=0.85, capsize=4)
for i, v in enumerate(res.auc_mean):
    ax.text(i, v + 0.02, f"{v:.3f}", ha="center", fontsize=10)
ax.set_ylim(0.7, 1.02); ax.set_ylabel("AUC (10×5-fold CV)")
ax.set_title("Diagnostic performance of four machine-learning models"); ax.grid(alpha=0.3, axis="y")
ax.tick_params(axis="x", labelsize=9)

ax = axes[1]
show = imp.head(10).iloc[::-1]
xs = np.arange(len(show))
ax.barh(xs - 0.25, show.RF / show.RF.max(), 0.22, label="Random forest", color="#4a7c96")
ax.barh(xs, show.GBM / show.GBM.max(), 0.22, label="GBM", color="#b08a3e")
ax.barh(xs + 0.25, show.GLM / show.GLM.max(), 0.22, label="GLM |coef|", color="#7c8f4a")
ax.set_yticks(xs); ax.set_yticklabels(show.index, fontsize=9)
ax.set_xlabel("Normalized importance")
ax.set_title("Consensus feature-importance ranking (top 10)")
ax.legend(fontsize=8, loc="lower right"); ax.grid(alpha=0.3, axis="x")
fig.tight_layout()
fig.savefig(f"{OUT}/EN_figureB5_model_comparison.png", dpi=150, bbox_inches="tight")
print("figureB5 saved")

# ---------- 图2: ROC (签名 + 单基因) ----------
fig, ax = plt.subplots(figsize=(6.4, 6))
ax.plot(fpr_c, tpr_c, color="#c0392b", lw=2.2,
        label=f"5-gene signature (CV AUC={cv_auc:.3f})")
ax.plot(fpr_t, tpr_t, color="#c0392b", lw=1, ls="--", alpha=0.5, label=f"Training fit AUC={auc(fpr_t, tpr_t):.3f}")
cmap = plt.cm.tab10
for i, (_, r) in enumerate(single.iterrows()):
    xg = expr.loc[r.gene, train_cols].values.reshape(-1, 1)
    m = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000)).fit(xg, y)
    pg = m.predict_proba(xg)[:, 1]
    f, t, _ = roc_curve(y, pg)
    ax.plot(f, t, color=cmap(i + 2), lw=1.1, alpha=0.8, label=f"{r.gene} (AUC={r.cv_auc:.3f})")
ax.plot([0, 1], [0, 1], "k--", lw=0.8)
ax.set_xlabel("1 \u2212 specificity"); ax.set_ylabel("Sensitivity")
ax.set_title("ROC of the senescence gene signature for i-IFTA")
ax.legend(fontsize=8.5, loc="lower right"); ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(f"{OUT}/EN_figureB6_ROC.png", dpi=150, bbox_inches="tight")
print("figureB6 saved")

# ---------- 图3: 三组评分谱 (IFTA 不入模, 作谱系验证) ----------
all_score = sig.predict_proba(expr.loc[HUB].T.values)[:, 1]
S = pd.DataFrame({"score": all_score, "group": grp.values})
fig, ax = plt.subplots(figsize=(7, 4.8))
order = ["Norm", "IFTA", "IFTA+I"]
disp = {"Norm": "Norm", "IFTA": "pure IF/TA", "IFTA+I": "i-IFTA"}
bp = ax.boxplot([S.loc[S.group == g, "score"] for g in order],
                tick_labels=[f"{disp[g]}\n(n={int((S.group==g).sum())})" for g in order],
                patch_artist=True, widths=0.5)
for patch, c in zip(bp["boxes"], ["#4a7c96", "#b08a3e", "#c96a4a"]):
    patch.set_facecolor(c); patch.set_alpha(0.75)
for i, g in enumerate(order):
    xs = np.random.RandomState(1).normal(i + 1, 0.05, size=int((S.group == g).sum()))
    ax.scatter(xs, S.loc[S.group == g, "score"], s=14, color="k", alpha=0.5, zorder=3)
p_ni = stats.mannwhitneyu(S.loc[S.group == "Norm", "score"], S.loc[S.group == "IFTA+I", "score"]).pvalue
p_nf = stats.mannwhitneyu(S.loc[S.group == "Norm", "score"], S.loc[S.group == "IFTA", "score"]).pvalue
p_fi = stats.mannwhitneyu(S.loc[S.group == "IFTA", "score"], S.loc[S.group == "IFTA+I", "score"]).pvalue
def fmt_p(p):
    return "P < 0.001" if p < 0.001 else f"P = {p:.3f}"
ax.set_title(f"Distribution of the 5-gene signature across all 65 biopsies\nNorm vs pure IF/TA: {fmt_p(p_nf)} | pure IF/TA vs i-IFTA: {fmt_p(p_fi)} | Norm vs i-IFTA: {fmt_p(p_ni)}")
ax.set_ylabel("Predicted probability of i-IFTA"); ax.grid(alpha=0.3, axis="y")
fig.tight_layout()
fig.savefig(f"{OUT}/EN_figureB7_signature_distribution.png", dpi=150, bbox_inches="tight")
print("figureB7 saved")

# ---------- 图4: 列线图 + 校准 + DCA (手绘, 基于GLM系数) ----------
lr2 = LogisticRegression(max_iter=2000).fit(Xh, y)   # 未标准化, 用于可解释列线图
beta = lr2.coef_[0]; b0 = lr2.intercept_[0]
mx = Xh.max(axis=0); mn = Xh.min(axis=0)
lpmin = (mn * np.where(beta >= 0, beta, 0)).sum() + (mx * np.where(beta < 0, beta, 0)).sum() + b0
lpmax = (mx * np.where(beta >= 0, beta, 0)).sum() + (mn * np.where(beta < 0, beta, 0)).sum() + b0
span = max(abs(beta * (mx - mn))); pts_per = beta * (mx - mn) / span * 100
T_max = float(np.abs(pts_per).sum())     # 理论最大总分
ax_right = min(T_max, 500)

fig, ax = plt.subplots(figsize=(9.5, 4.6))
n = len(HUB)
ax.set_xlim(-8, ax_right + 8); ax.set_ylim(-0.8, n + 2.6); ax.axis("off")
# 分值轴
for frac in [0, 0.25, 0.5, 0.75, 1]:
    xv = frac * 100
    ax.plot([xv, xv], [-0.3, n - 0.6], color="#dddddd", lw=0.6, zorder=0)
    ax.text(xv, -0.55, f"{xv:.0f}", ha="center", fontsize=8, color="#555555")
ax.text(50, -0.15, "分值", ha="center", fontsize=9, color="#555555")
for i, g in enumerate(HUB):
    yv = n - 1 - i
    ax.text(-6, yv, g, ha="right", va="center", fontsize=10, fontweight="bold")
    lo, hi = (mn[i], mx[i]) if beta[i] >= 0 else (mx[i], mn[i])
    ax.plot([0, 100], [yv, yv], color="#4a7c96", lw=1.6)
    ax.scatter([0, 100], [yv, yv], s=18, color="#c0392b", zorder=3)
    ax.text(0, yv + 0.22, f"{lo:.1f}", ha="center", fontsize=7.5, color="#c0392b")
    ax.text(100, yv + 0.22, f"{hi:.1f}", ha="center", fontsize=7.5, color="#c0392b")
ax.text(-6, n - 0.5, "总分", ha="right", va="center", fontsize=10, fontweight="bold")
ax.plot([0, ax_right], [n - 0.5, n - 0.5], color="#b08a3e", lw=2.2)
for frac in [0, 0.25, 0.5, 0.75, 1]:
    xv = frac * ax_right
    ax.scatter([xv], [n - 0.5], s=12, color="#8a6d2f", zorder=3)
    ax.text(xv, n - 0.18, f"{xv:.0f}", ha="center", fontsize=7.5, color="#8a6d2f")
# 风险轴: lp = lpmin + 总分/100*span → 按轴比例反算概率
ax.text(-6, n + 0.7, "风险", ha="right", va="center", fontsize=10, fontweight="bold")
ax.plot([0, ax_right], [n + 0.7, n + 0.7], color="#c96a4a", lw=2.2)
for frac in [0, 0.25, 0.5, 0.75, 1]:
    xv = frac * ax_right
    lp = lpmin + xv / 100 * span
    pv = 1 / (1 + np.exp(-lp))
    ax.scatter([xv], [n + 0.7], s=16, color="#8A5A44", zorder=3)
    ax.text(xv, n + 1.05, f"{pv:.2f}", ha="center", fontsize=7.5, color="#8A5A44")
ax.set_title("Nomogram of the 5-gene signature for IF/TA+I", fontsize=12)
fig.tight_layout()
fig.savefig(f"{OUT}/EN_figureB8_nomogram.png", dpi=150, bbox_inches="tight")
print("figureB8 saved")

# 校准曲线 (分组均值 vs 预测)
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
ax = axes[0]
nb = 4
bins = np.quantile(p_all, np.linspace(0, 1, nb + 1)); bins[0], bins[-1] = -0.01, 1.01
cal = pd.DataFrame({"p": p_all, "y": y})
cal["bin"] = np.digitize(cal.p, bins) - 1
gc = cal.groupby("bin").agg(pred=("p", "mean"), obs=("y", "mean"), n=("y", "size"))
ax.plot([0, 1], [0, 1], "k--", lw=0.8, label="Ideal")
ax.plot(gc.pred, gc.obs, "o-", color="#c0392b", lw=1.8, label="Observed")
for _, r in gc.iterrows():
    ax.annotate(f"n={int(r.n)}", (r.pred, r.obs), textcoords="offset points",
                xytext=(6, 6), fontsize=8)
ax.set_xlabel("Predicted probability"); ax.set_ylabel("Observed frequency")
ax.set_title(f"Calibration curve (CV AUC={cv_auc:.3f})"); ax.legend(fontsize=9); ax.grid(alpha=0.3)

ax = axes[1]
ths = np.linspace(0.02, 0.98, 60)
n_pos, n_neg = y.sum(), (1 - y).sum()
net_benefit_model, net_benefit_all = [], []
for t in ths:
    pred_pos = p_all >= t
    tp = (pred_pos & (y == 1)).sum(); fp = (pred_pos & (y == 0)).sum()
    w = t / (1 - t)
    net_benefit_model.append(tp / n_pos - fp / n_neg * w)
    net_benefit_all.append(1 - (n_neg / n_pos) * w)
ax.plot(ths, net_benefit_model, color="#c0392b", lw=1.8, label="5-gene signature")
ax.plot(ths, net_benefit_all, color="#777777", lw=1.2, ls="--", label="Treat all")
ax.axhline(0, color="#4a7c96", lw=1.2, ls="-", label="Treat none")
ax.set_ylim(-0.3, 1.0); ax.set_xlabel("Threshold probability"); ax.set_ylabel("Net benefit")
ax.set_title("Decision curve analysis (DCA)"); ax.legend(fontsize=9); ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(f"{OUT}/EN_figureB9_calibration_DCA.png", dpi=150, bbox_inches="tight")
print("figureB9 saved")

# 保存签名细节
pd.DataFrame({"gene": HUB, "beta": beta, "points_range": mx - mn}).to_csv(
    f"{OUT}/40_EN_signature_coefficients.csv", index=False, encoding="utf-8-sig")
S.to_csv(f"{OUT}/40_EN_sample_scores.csv", encoding="utf-8-sig")
print("done")
