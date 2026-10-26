# -*- coding: utf-8 -*-
"""
论文B线 脚本41:
  Part A: 人肾图谱2.5 (304,989细胞核) 中 21个ARDEGs 的细胞类型定位
  Part B: GSE22459 中 5基因签名的免疫混杂解析 (与区室标志相关 + 偏相关)
用法:
  python 41_图谱定位与免疫解析.py chunks <start> <end>   # 流式读图谱表达
  python 41_图谱定位与免疫解析.py analyze                # 统计+作图
"""
import sys, os, time
import numpy as np
import h5py
from scipy import sparse

H5 = "geo_data/02d5db96-d67e-45f2-b556-b6bd99753665.h5ad"
STATE = "results/41_state.npz"
OUT = "results/fibrosis_ARG"
NCHUNK = 40000
N = 304989

import pandas as pd
GENES = pd.read_csv(f"{OUT}/39_ARDEGs.csv", index_col=0).index.tolist()
print("ARDEGs:", len(GENES))

def gene_indices(f):
    fn = f["var"]["feature_name"]
    cats = fn["categories"][:].astype(str); codes = fn["codes"][:]
    idx = []
    for g in GENES:
        hit = np.where(cats == g)[0]
        idx.append(np.where(codes == hit[0])[0][0] if len(hit) else -1)
    return np.array(idx)

def run_chunks(a, b):
    t0 = time.time()
    f = h5py.File(H5, "r")
    X = f["X"]; indptr = X["indptr"]
    if os.path.exists(STATE):
        st = dict(np.load(STATE))
        vals = st["vals"].copy(); sums = st["sums"].copy(); done = set(int(x) for x in st["done"])
    else:
        vals = np.full((N, len(GENES)), np.nan, dtype=np.float32)
        sums = np.zeros(len(GENES)); done = set()
    gidx = gene_indices(f)
    starts = list(range(0, N, NCHUNK))
    for ci in range(a, min(b, len(starts))):
        if ci in done: continue
        if time.time() - t0 > 230: print("时间到,保存退出"); break
        r0, r1 = starts[ci], min(starts[ci] + NCHUNK, N)
        p0, p1 = int(indptr[r0]), int(indptr[r1])
        data = X["data"][p0:p1]; indices = X["indices"][p0:p1]
        ip = indptr[r0:r1 + 1].astype(np.int64) - p0
        chunk = sparse.csr_matrix((data, indices, ip), shape=(r1 - r0, 35439)).tocsc()
        sub = np.asarray(chunk[:, gidx[gidx >= 0]].todense())
        vals[r0:r1][:, gidx >= 0] = sub
        sums += sub.sum(axis=0)
        done.add(ci)
        print(f"chunk {ci} done {r0}-{r1} ({time.time()-t0:.0f}s)", flush=True)
        np.savez(STATE, vals=vals, sums=sums, done=np.array(sorted(done)), gidx=gidx)
        del chunk, sub, data, indices
    print(f"progress {len(done)}/{len(starts)}")

def run_analyze():
    from scipy import stats
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from pathlib import Path
    sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
    from daimon_runtime import setup_plot
    setup_plot()
    os.makedirs(OUT, exist_ok=True)

    f = h5py.File(H5, "r")
    obs = f["obs"]
    def cat(k):
        d = obs[k]; return np.asarray(d["categories"][:].astype(str))[d["codes"][:]]
    l1 = cat("subclass.l1")
    donor = cat("donor_id")
    dev = cat("development_stage")

    st = np.load(STATE)
    vals = st["vals"]; gidx = st["gidx"]
    ok = gidx >= 0
    print(f"图谱中匹配到 {ok.sum()}/{len(GENES)} 基因")
    g_present = [g for g, o in zip(GENES, ok) if o]
    V = vals[:, ok]                     # (N, g)

    # ---- Part A: 细胞类型定位 (dot plot: 表达率 + 平均表达) ----
    L1_ZH = {"IMM": "Immune", "FIB": "Fibroblast", "VSM/P": "Vascular SM/pericyte",
             "EC": "Endothelial", "PT": "Proximal tubule", "POD": "Podocyte", "PEC": "Parietal epi.",
             "DTL": "Descending thin limb", "ATL": "Ascending thin limb", "DCT": "Distal tubule", "TAL": "Thick ascend. limb",
             "IC": "Intercalated", "CNT": "Connecting tubule", "PC": "Principal cell",
             "PapE": "Papillary epi.", "NEU": "Neutrophil"}
    ct_order = sorted(pd.Series(l1).value_counts().index.tolist(), key=lambda c: -pd.Series(l1).value_counts()[c])
    frac = np.zeros((len(g_present), len(ct_order)))
    mean_expr = np.zeros_like(frac)
    for j, ct in enumerate(ct_order):
        m = l1 == ct
        sub = V[m]
        frac[:, j] = (sub > 0).mean(axis=0)
        mean_expr[:, j] = np.nanmean(np.where(sub > 0, sub, np.nan), axis=0)
    F = pd.DataFrame(frac, index=g_present, columns=ct_order)
    M = pd.DataFrame(mean_expr, index=g_present, columns=ct_order)
    F.to_csv(f"{OUT}/41_EN_cell_fraction.csv", encoding="utf-8-sig")
    M.to_csv(f"{OUT}/41_EN_cell_meanexpr.csv", encoding="utf-8-sig")

    # 每个基因的主表达细胞类型
    top_ct = F.idxmax(axis=1)
    print("\n每个ARDEG表达最高的细胞类型:")
    for g in g_present:
        print(f"  {g}: {top_ct[g]} (表达率 {F.loc[g, top_ct[g]]*100:.0f}%)")
    pd.DataFrame({"gene": g_present, "top_celltype": top_ct,
                  "frac": [F.loc[g, top_ct[g]] for g in g_present]}).to_csv(
        f"{OUT}/41_EN_top_celltype.csv", index=False, encoding="utf-8-sig")

    # dot plot
    fig, ax = plt.subplots(figsize=(0.55 * len(ct_order) + 3, 0.32 * len(g_present) + 2))
    xs = np.arange(len(ct_order)); ys = np.arange(len(g_present))
    for i, g in enumerate(g_present):
        for j, ct in enumerate(ct_order):
            fv = F.loc[g, ct]; mv = M.loc[g, ct]
            if fv <= 0.02: continue
            ax.scatter(j, len(g_present) - 1 - i, s=30 + fv * 420,
                       c=[mv], cmap="viridis", vmin=0, vmax=np.nanmax(M.values), alpha=0.9)
    ax.set_xticks(xs); ax.set_xticklabels([L1_ZH.get(c, c) for c in ct_order], rotation=45, ha="right", fontsize=9)
    ax.set_yticks(ys[::-1]); ax.set_yticklabels(g_present, fontsize=9)
    ax.set_title("Expression of 21 ARDEGs across cell types in the human kidney atlas v1.5\ndot size = fraction of expressing cells, colour = mean expression")
    fig.tight_layout()
    fig.savefig(f"{OUT}/EN_figureB10_cell_localization.png", dpi=150, bbox_inches="tight")
    print("figureB10 saved")

    # ---- Part B: 免疫混杂解析 (GSE22459 bulk) ----
    expr = pd.read_csv(f"{OUT}/39_GSE22459_log2expr.csv", index_col=0)
    grp = pd.read_csv(f"{OUT}/39_GSE22459_groups.csv", index_col=0).iloc[:, 0]
    sig = pd.read_csv(f"{OUT}/40_all_sample_scores.csv", index_col=0)
    z = expr.sub(expr.mean(axis=1), axis=0).div(expr.std(axis=1).replace(0, np.nan), axis=0)
    COMP = {"T cells": ["CD3D", "CD3E", "CD2", "TRAC"],
            "Macrophages": ["CD68", "CD14", "LST1", "FCGR3A"],
            "B cells": ["CD79A", "MS4A1"],
            "NK cells": ["NKG7", "GNLY"],
            "Fibroblast/Stroma": ["DCN", "LUM", "PDGFRB", "COL1A2"],
            "Proximal tubule": ["LRP2", "CUBN", "SLC34A1"]}
    comp_scores = {}
    for k, gs in COMP.items():
        avail = [g for g in gs if g in z.index]
        comp_scores[k] = z.loc[avail].mean(axis=0)
    C = pd.DataFrame(comp_scores)
    C["signature"] = sig["score"].values        # 位置对齐 (签名CSV为RangeIndex)
    C["group"] = grp.values
    C.to_csv(f"{OUT}/41_EN_signature_compartments.csv", encoding="utf-8-sig")

    def pcorr(a, b, ctrl):
        D = pd.concat([a, b] + ctrl, axis=1).dropna()
        R = D.corr().values
        if np.isnan(R).any(): return np.nan, len(D)
        Ci = np.linalg.pinv(R)
        return -Ci[0, 1] / np.sqrt(Ci[0, 0] * Ci[1, 1]), len(D)

    rank = grp.map({"Norm": 0, "IFTA": 1, "IFTA+I": 2}).astype(float)
    res_rows = []
    for k in COMP:
        r_raw = C[k].corr(C["signature"])
        r_adj, n = pcorr(C["signature"].rename("s"), C[k].rename("c"), [])
        res_rows.append((k, r_raw))
        print(f"签名 vs {k}: r = {r_raw:.2f}")
    # 偏相关: 签名~IF/TA等级 | 控制T细胞+巨噬
    r_crude, _ = pcorr(C["signature"].rename("s"), rank.rename("r"), [])
    r_adj, _ = pcorr(C["signature"].rename("s"), rank.rename("r"),
                     [C["T cells"].rename("t"), C["Macrophages"].rename("m")])
    print(f"\n[来源纪律] 签名 vs IF/TA等级: r原始={r_crude:.2f} → 控制T细胞+巨噬后 r={r_adj:.2f}")

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6))
    ax = axes[0]
    ks = list(COMP.keys())
    rs = [C[k].corr(C["signature"]) for k in ks]
    bars = ax.bar(ks, rs, color=["#c96a4a", "#c96a4a", "#c96a4a", "#c96a4a", "#b08a3e", "#4a7c96"], alpha=0.85)
    for b_, v in zip(bars, rs):
        ax.text(b_.get_x() + b_.get_width() / 2, v + 0.02, f"{v:.2f}", ha="center", fontsize=9)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_ylabel("Correlation r with the 5-gene signature")
    ax.set_title("Correlation of the 5-gene signature with compartment markers\n(immune compartments → high; tubular epithelium → low)")
    ax.tick_params(axis="x", labelsize=8, rotation=20); ax.grid(alpha=0.3, axis="y")

    ax = axes[1]
    labels = ["Crude", "Adj. T cells + macro."]
    vals_ = [r_crude, r_adj]
    bars = ax.bar(labels, vals_, color=["#4a7c96", "#c96a4a"], alpha=0.85, width=0.45)
    for b_, v in zip(bars, vals_):
        ax.text(b_.get_x() + b_.get_width() / 2, v + 0.02, f"{v:.2f}", ha="center", fontsize=11)
    ax.set_ylabel("Partial correlation r with IF/TA grade")
    ax.set_title(f"Source-discipline test: correlation with lesion grade\nretained after adjusting for immune infiltrates (r={r_adj:.2f})")
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(f"{OUT}/EN_figureB11_immune_deconvolution.png", dpi=150, bbox_inches="tight")
    print("figureB11 saved")
    print("done")

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "chunks":
        run_chunks(int(sys.argv[2]), int(sys.argv[3]))
    else:
        run_analyze()
