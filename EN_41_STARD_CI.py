# -*- coding: utf-8 -*-
"""
脚本50 (STARD items 23/24): 主签名5折CV AUC的bootstrap 95%CI
复刻 EN_40.py 的完整建模流程(random_state固定,可复现),验证锁定数字0.907/0.918,
再对41个样本的out-of-fold预测做2000次bootstrap → 百分位95%CI
另附: 队列中位数切点下的敏感度/特异度(item 23 探索性交叉表)
"""
import numpy as np, pandas as pd
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score

OUT = "results/纤维化ARG"
expr = pd.read_csv(f"{OUT}/39_GSE22459_log2expr.csv", index_col=0)
grp = pd.read_csv(f"{OUT}/39_GSE22459_groups.csv", index_col=0).iloc[:, 0]
genes = pd.read_csv(f"{OUT}/39_ARDEGs.csv", index_col=0).index.tolist()

train_mask = grp.isin(["Norm", "IFTA+I"])
train_cols = grp[train_mask].index
X = expr.loc[genes, train_cols].T.values
y = (grp[train_mask] == "IFTA+I").astype(int).values
names = np.array(genes)

# ---- 复刻特征重要性共识(与EN_40完全一致) ----
imp = pd.DataFrame(index=names)
rf = RandomForestClassifier(n_estimators=1000, random_state=42).fit(X, y)
imp["RF"] = rf.feature_importances_
lr = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000)).fit(X, y)
imp["GLM"] = np.abs(lr.named_steps["logisticregression"].coef_[0])
gb = GradientBoostingClassifier(random_state=42).fit(X, y)
imp["GBM"] = gb.feature_importances_
for c in ["RF", "GLM", "GBM"]:
    imp[c + "_rank"] = imp[c].rank(ascending=False)
imp["mean_rank"] = imp[[c + "_rank" for c in ["RF", "GLM", "GBM"]]].mean(axis=1)
imp = imp.sort_values("mean_rank")
HUB = imp.head(5).index.tolist()
print("HUB(top5):", HUB)
assert HUB == ["TLR10", "RRM2", "TOP2A", "TLR2", "GRPR"], f"HUB与锁定不符: {HUB}"

cv = StratifiedKFold(5, shuffle=True, random_state=42)

def sig_cv_auc(hub):
    Xh = expr.loc[hub, train_cols].T.values
    sig = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000)).fit(Xh, y)
    cv_p = cross_val_predict(sig, Xh, y, cv=cv, method="predict_proba")[:, 1]
    return roc_auc_score(y, cv_p), cv_p

def boot_ci(y_true, p, n_boot=2000, seed=42):
    rng = np.random.RandomState(seed)
    n = len(y_true)
    aucs = []
    for _ in range(n_boot):
        idx = rng.randint(0, n, n)
        if len(np.unique(y_true[idx])) < 2:
            continue
        aucs.append(roc_auc_score(y_true[idx], p[idx]))
    aucs = np.array(aucs)
    return np.percentile(aucs, 2.5), np.percentile(aucs, 97.5)

rows = []
oof_all = {}
for label, hub in [("5-gene signature", HUB), ("4-gene signature (TLR10-free)", HUB[1:])]:
    auc_cv, cv_p = sig_cv_auc(hub)
    lo, hi = boot_ci(y, cv_p)
    # 探索性交叉表: 队列中位数切点
    cut = np.median(cv_p)
    pred_pos = cv_p >= cut
    tp = int((pred_pos & (y == 1)).sum()); fp = int((pred_pos & (y == 0)).sum())
    fn = int((~pred_pos & (y == 1)).sum()); tn = int((~pred_pos & (y == 0)).sum())
    sens, spec = tp / (tp + fn), tn / (tn + fp)
    print(f"{label}: CV AUC={auc_cv:.3f}, boot 95%CI=({lo:.3f}-{hi:.3f}), "
          f"median-cut sens={sens:.2f} spec={spec:.2f} (TP{tp}/FP{fp}/FN{fn}/TN{tn})")
    rows.append([label, f"+".join(hub), round(auc_cv, 3), round(lo, 3), round(hi, 3),
                 tp, fp, fn, tn, round(sens, 3), round(spec, 3)])
    oof_all[label] = (cv_p, y)

pd.DataFrame(rows, columns=["signature", "genes", "cv_auc_5fold", "boot95_ci_low", "boot95_ci_high",
                            "TP", "FP", "FN", "TN", "sensitivity_median_cut", "specificity_median_cut"]
             ).to_csv(f"{OUT}/41_EN_CV_AUC_bootstrap_CI.csv", index=False, encoding="utf-8-sig")

oof = pd.DataFrame({"y_true": y, "oof_5gene": oof_all["5-gene signature"][0],
                    "oof_4gene": oof_all["4-gene signature (TLR10-free)"][0]}, index=train_cols)
oof.to_csv(f"{OUT}/41_EN_oof_predictions.csv", encoding="utf-8-sig")
print("saved: 41_EN_CV_AUC_bootstrap_CI.csv, 41_EN_oof_predictions.csv")
