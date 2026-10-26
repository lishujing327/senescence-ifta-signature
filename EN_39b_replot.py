# -*- coding: utf-8 -*-
"""39b 英文重绘：PPI 网络 + 富集气泡图。
离线版：直接读已缓存的 39b_PPI_edges.csv / 39b_STRING_enrichment_full.csv / 39_ARDEGs.csv，
不调用 STRING API。图形逻辑与 39b_ARG_PPI富集.py 完全一致（圆形布局、score≥0.4、各取top8）。
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
from daimon_runtime import setup_plot
setup_plot()

OUT = "results/fibrosis_ARG"

net = pd.read_csv(f"{OUT}/39b_PPI_edges.csv")
enr = pd.read_csv(f"{OUT}/39b_STRING_enrichment_full.csv")
deg = pd.read_csv(f"{OUT}/39_ARDEGs.csv", index_col=0)
print(f"边表 {len(net)} 条（缓存已按 score≥0.4 过滤）；富集 {len(enr)} 条")

# ---------- PPI 网络图 ----------
nodes = sorted(set(net.preferredName_A) | set(net.preferredName_B))
pos = {n: (np.cos(2 * np.pi * i / len(nodes)), np.sin(2 * np.pi * i / len(nodes)))
       for i, n in enumerate(nodes)}

fig, ax = plt.subplots(figsize=(8.2, 8.2))
for _, r in net.iterrows():
    a, b = r.preferredName_A, r.preferredName_B
    ax.plot([pos[a][0], pos[b][0]], [pos[a][1], pos[b][1]],
            color="#999999", lw=0.4 + r.score * 3, alpha=0.7, zorder=1)
for n in nodes:
    x, y = pos[n]
    fc = deg.loc[n, "log2FC"] if n in deg.index else 0
    color = "#c0392b" if fc > 0 else "#2c5f8a"
    ax.scatter(x, y, s=600 + 400 * abs(fc), c=color, alpha=0.85, zorder=2,
               edgecolors="white", linewidths=1.2)
    ax.annotate(n, (x, y), textcoords="offset points", xytext=(0, 11),
                ha="center", fontsize=9, fontweight="bold")
ax.set_title("PPI network of ARDEGs (STRING, score ≥ 0.4)\nred = upregulated, blue = downregulated, node size = |log2FC|", fontsize=12)
ax.axis("off")
fig.tight_layout()
fig.savefig(f"{OUT}/EN_figureB3_PPI_network.png", dpi=150, bbox_inches="tight")
print("EN_figureB3 saved")

# ---------- 富集气泡图 ----------
picks = []
for cat in ["GO Biological Process", "KEGG", "Reactome"]:
    sub = enr[enr.category == cat].sort_values("fdr").head(8)
    for _, r in sub.iterrows():
        picks.append((cat, r["term"], r["fdr"], r["number_of_genes"]))
B = pd.DataFrame(picks, columns=["cat", "term", "fdr", "n"]).sort_values("fdr", ascending=False)
fig, ax = plt.subplots(figsize=(10, 0.42 * len(B) + 1.6))
ypos = np.arange(len(B))
sc = ax.scatter(-np.log10(B.fdr.clip(lower=1e-300)), ypos, s=40 + B.n * 22,
                c=B["n"], cmap="viridis", alpha=0.85, edgecolors="k", linewidths=0.4)
ax.set_yticks(ypos)
ax.set_yticklabels([f"[{c.split()[-1]}] {t[:62]}" for c, t in zip(B.cat, B.term)], fontsize=8.5)
ax.set_xlabel("-log10(FDR)")
ax.set_title("Functional enrichment of ARDEGs (STRING)")
fig.colorbar(sc, ax=ax, shrink=0.7, label="Gene count")
ax.grid(alpha=0.25, axis="x")
fig.tight_layout()
fig.savefig(f"{OUT}/EN_figureB4_enrichment.png", dpi=150, bbox_inches="tight")
print("EN_figureB4 saved")
