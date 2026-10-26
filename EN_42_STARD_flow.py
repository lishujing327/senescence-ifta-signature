# -*- coding: utf-8 -*-
"""
STARD item 19: participant/dataset flow diagram (English, publication style)
Numbers all locked from manuscript EN_v2 / Supplementary Table S1
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from daimon_runtime import setup_plot
setup_plot()

C_MAIN = "#4a7c96"; C_VAL = "#7c8f4a"; C_OUT = "#8A5A44"; C_REF = "#b08a3e"

fig, ax = plt.subplots(figsize=(9.2, 10.5))
ax.set_xlim(0, 10); ax.set_ylim(0, 13); ax.axis("off")

def box(x, y, w, h, text, fc, ec=None, fs=9.5, bold=False, tc="white"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.12",
                                fc=fc, ec=ec or fc, lw=1.2))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
            color=tc, fontweight="bold" if bold else "normal", linespacing=1.5)

def arrow(x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                 mutation_scale=16, color="#555555", lw=1.4))

# Top source
box(1.8, 12.0, 6.4, 0.85, "Public transcriptomic datasets (GEO; CZ CELLxGENE)", "#35586b", fs=10.5, bold=True)

# GSE22459 branch
box(0.4, 10.3, 4.4, 1.2, "GSE22459 (GPL570)\n65 low-risk living-donor\n1-year protocol biopsies", C_MAIN)
box(5.2, 10.3, 4.4, 1.2, "Histological groups\nNorm 25  |  pure IF/TA 24  |  i-IFTA 16", C_MAIN, fs=9)
arrow(4.8, 10.9, 5.2, 10.9)

box(0.4, 8.4, 4.4, 1.1, "Model building (n = 41)\nNorm (25) + i-IFTA (16)\n10x5-fold CV, AUC metric", C_MAIN)
box(5.2, 8.4, 4.4, 1.1, "Lineage-gradient verification\npure IF/TA (n = 24)\nreserved, not modelled", C_REF)
arrow(2.6, 10.3, 2.6, 9.55)
arrow(7.4, 10.3, 7.4, 9.55)

# Feature selection + signature
box(0.4, 6.5, 4.4, 1.2, "Consensus feature importance\n(RF + GLM + GBM)\ntop-5 hub genes", C_MAIN)
box(5.2, 6.5, 4.4, 1.2, "5-gene signature\nTLR10, RRM2, TOP2A, TLR2, GRPR\n4-gene portable version", C_MAIN)
arrow(2.6, 8.4, 2.6, 7.75)
arrow(4.8, 7.1, 5.2, 7.1)

# Single-nucleus atlas
box(0.4, 4.6, 4.4, 1.1, "Human Kidney snRNA-seq\nAtlas v1.5 (304,989 nuclei)\n-> single-nucleus localization", "#6b4a7c")
arrow(2.6, 6.5, 2.6, 5.75)

# External validation
box(5.2, 4.6, 4.4, 1.1, "External validation (no refitting)\nGSE53605 (n = 55, GPL571)\nGSE7392 (30 paired, GPL570)", C_VAL)
arrow(7.4, 6.5, 7.4, 5.75)

# Prognostic cohort
box(0.4, 2.7, 9.2, 1.1, "Prognostic testing (no refitting)\nGSE21374 (GPL570): 105 late for-cause biopsies, 30 graft-loss events\nCox regression per 1-SD increase + median-split Kaplan-Meier", C_VAL, fs=9.5)
arrow(2.6, 4.6, 2.6, 3.85)
arrow(7.4, 4.6, 7.4, 3.85)

# Outputs
box(0.4, 0.7, 9.2, 1.3, "Locked results\n5-gene signature 5-fold CV AUC 0.907 (95% CI 0.784-0.989); 4-gene 0.918\nGSE53605 IF/TA AUC 1.000; GSE21374 graft loss HR 1.57 (1.19-2.07)\nBenchmarking vs whole-genome/random sets; overlap with 5 classifiers = 0", C_OUT, fs=9.5)
arrow(5.0, 2.7, 5.0, 2.05)

ax.set_title("Flow of datasets and analyses (STARD item 19)", fontsize=12, fontweight="bold", pad=14)
fig.savefig("results/纤维化ARG/EN_figureB16_STARD_flow.png", dpi=300, bbox_inches="tight")
print("saved: EN_figureB16_STARD_flow.png")
