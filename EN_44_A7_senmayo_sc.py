# -*- coding: utf-8 -*-
"""
A7: SenMayo senescence scoring on the human kidney single-nucleus atlas
(purely computational, zero cost).
- SenMayo gene set: Saul et al. Nat Commun 2022;13:4827 (125 genes;
  deliberately excludes CDKN2A/CDKN1A).
- Atlas: CZ CELLxGENE human kidney atlas v1.5 (304,989 nuclei, 35,439 genes,
  CSR sparse h5ad).
- Question: which cell type carries the senescence load, and does PT-cell
  senescence increase with donor age?
Chunked CSR reading controls memory; the score is the per-cell mean of
cross-cell z values of target-gene expression (simplified AddModuleScore).
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
import h5py, numpy as np, pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from daimon_runtime import setup_plot
setup_plot()

OUT = "results/fibrosis_ARG"
H5 = "geo_data/02d5db96-d67e-45f2-b556-b6bd99753665.h5ad"

SENMAYO = """ACVR1B CCL7 EGFR IGFBP6 MMP1 SCAMP4 ANG CCL8 EREG IGFBP7 MMP10 SELPLG
ANGPT1 CD55 ESM1 IL10 MMP12 SEMA3F ANGPTL4 CD9 ETS2 IL13 MMP13 SERPINB4
AREG CSF1 FAS IL15 MMP14 SERPINE1 AXL CSF2 FGF1 IL18 MMP2 SERPINE2
BEX3 CSF2RB FGF2 IL1A MMP3 SPP1 BMP2 CST4 FGF7 IL1B MMP9 SPX
BMP6 CTNNB1 GDF15 IL2 NAP1L4 TIMP2 C3 CTSB GEM IL32 NRG1 TNF
CCL1 CXCL1 GMFG IL6 PAPPA TNFRSF10C CCL13 CXCL10 HGF IL6ST PECAM1 TNFRSF11B
CCL16 CXCL12 HMGB1 IL7 PGF TNFRSF1A CCL2 CXCL16 ICAM1 INHA PIGF TNFRSF1B
CCL20 CXCL2 ICAM3 IQGAP2 PLAT TUBGCP2 CCL24 CXCL3 IGF1 ITGA2 PLAU VEGFA
CCL26 CXCL8 IGFBP1 ITPKA PLAUR VEGFC CCL3 CXCR2 IGFBP2 JUN PTBP1 VGF
CCL3L1 DKK1 IGFBP3 KITLG PTGER2 WNT16 CCL4 EDN1 IGFBP4 LCP1 PTGES WNT2
CCL5 EGF IGFBP5 MIF RPS6KA5""".split()
print("SenMayo genes:", len(SENMAYO))

f = h5py.File(H5, "r")
X = f["X"]
indptr = X["indptr"][:]
n_rows, n_cols = X.attrs["shape"]
print("X:", n_rows, "x", n_cols, "| indptr ok")

def cat_to_array(node):
    cats = np.array([c.decode() if isinstance(c, bytes) else c for c in node["categories"][:]])
    codes = node["codes"][:]
    return pd.Categorical.from_codes(codes, categories=cats)

def str_to_array(node):
    # var feature_name may be categorical or variable-length string
    if isinstance(node, h5py.Group):
        return np.array(cat_to_array(node))
    return np.array([x.decode() if isinstance(x, bytes) else x for x in node[:]])

genes = str_to_array(f["var"]["feature_name"])
obs = pd.DataFrame({
    "cell_type": cat_to_array(f["obs"]["cell_type"]),
    "subclass": cat_to_array(f["obs"]["subclass.l1"]),
    "disease_cat": cat_to_array(f["obs"]["disease_category"]),
    "dev_stage": cat_to_array(f["obs"]["development_stage"]),
    "donor": cat_to_array(f["obs"]["donor_id"]),
})
targets = sorted(set(SENMAYO) | {"CDKN2A", "CDKN1A"})
col_idx = {g: np.where(genes == g)[0] for g in targets}
found = {g: int(v[0]) for g, v in col_idx.items() if len(v)}
print(f"target genes found on platform: {len(found)}/127 (SenMayo: {len([g for g in found if g in SENMAYO])}/125)")
sen_present = [g for g in SENMAYO if g in found]
col_list = np.array([found[g] for g in sen_present + ["CDKN2A", "CDKN1A"] if g in found])
gnames = [g for g in sen_present + ["CDKN2A", "CDKN1A"] if g in found]

# ---------- chunked reading of target columns ----------
CHUNK = 20000
data_by_gene = {g: np.zeros(n_rows, dtype=np.float32) for g in gnames}
pos = 0
for s in range(0, n_rows, CHUNK):
    e = min(s + CHUNK, n_rows)
    d = X["data"][indptr[s]:indptr[e]]
    ix = X["indices"][indptr[s]:indptr[e]]
    rows_local = np.repeat(np.arange(e - s), np.diff(indptr[s:e + 1]))
    for j, gname in enumerate(gnames):
        c = col_list[j]
        m = ix == c
        data_by_gene[gname][rows_local[m] + s] = d[m]
    if s % 100000 == 0:
        print(f"  rows {s}/{n_rows}")
X = None
f.close()

expr = pd.DataFrame(data_by_gene)
del data_by_gene
logexpr = np.log1p(expr)

# ---------- z scoring (simplified AddModuleScore: cross-cell z per target gene, mean per cell) ----------
z = (logexpr - logexpr.mean()) / logexpr.std(ddof=0)
score = z[sen_present].mean(axis=1)
obs["senmayo"] = score.values
obs["p16"] = logexpr["CDKN2A"].values if "CDKN2A" in logexpr else np.nan
obs["p21"] = logexpr["CDKN1A"].values if "CDKN1A" in logexpr else np.nan

# ---------- age mapping ----------
def stage_to_age(s):
    import re
    m = re.match(r"(\d+)-year-old", s)
    if m: return int(m.group(1))
    dec = {"second": 15, "third": 25, "fourth": 35, "fifth": 45, "sixth": 55, "seventh": 65, "eighth": 75}
    for k, v in dec.items():
        if k in s: return v
    return np.nan
obs["age"] = [stage_to_age(s) for s in obs["dev_stage"]]
print("age mapped:", obs["age"].notna().mean().round(3), "| range:", obs["age"].min(), "-", obs["age"].max())

# ---------- 1) SenMayo score by cell type ----------
ct = obs.groupby("cell_type", observed=True).agg(n=("senmayo", "size"), senmayo_mean=("senmayo", "mean"),
                                                 senmayo_sd=("senmayo", "std"),
                                                 p16_detect=("p16", lambda s: (s > 0).mean()),
                                                 p21_detect=("p21", lambda s: (s > 0).mean())).sort_values("senmayo_mean", ascending=False)
ct.to_csv(f"{OUT}/44_senmayo_by_celltype.csv", encoding="utf-8-sig")
print(ct.round(3).to_string())

# ---------- 2) PT cells: donor-age correlation ----------
pt = obs[obs["subclass"] == "PT"].copy()
healthy = pt[pt["disease_cat"].isin(["Healthy_reference_tissue", "Healthy_stone_donor"])]
donor_lvl = healthy.groupby("donor", observed=True).agg(age=("age", "median"), senmayo=("senmayo", "mean"),
                                                        n_cells=("senmayo", "size"),
                                                        disease_cat=("disease_cat", "first")).query("n_cells >= 20 & age == age")
r_h, p_h = stats.spearmanr(donor_lvl["age"], donor_lvl["senmayo"])
donor_lvl.to_csv(f"{OUT}/44_PT_senmayo_by_donor_age.csv", encoding="utf-8-sig")
print(f"\nHEALTHY PT donors (n={len(donor_lvl)}): SenMayo vs age Spearman r={r_h:.3f}, P={p_h:.4g}")

donor_all = pt.groupby("donor", observed=True).agg(age=("age", "median"), senmayo=("senmayo", "mean"),
                                                   n_cells=("senmayo", "size")).query("n_cells >= 20 & age == age")
r_a, p_a = stats.spearmanr(donor_all["age"], donor_all["senmayo"])
print(f"ALL PT donors (n={len(donor_all)}): r={r_a:.3f}, P={p_a:.4g}")

# reference: donor-level mean over ALL cells (not PT-specific)
donor_allcells = obs.groupby("donor", observed=True).agg(age=("age", "median"), senmayo=("senmayo", "mean")).query("age == age")
r_all, p_all = stats.spearmanr(donor_allcells["age"], donor_allcells["senmayo"])
print(f"ALL cells donor-level (n={len(donor_allcells)}): r={r_all:.3f}, P={p_all:.4g}")

# ---------- 3) PT cells: disease-state comparison ----------
dis = pt.groupby("disease_cat", observed=True).agg(n=("senmayo", "size"), senmayo_mean=("senmayo", "mean"), senmayo_sd=("senmayo", "std"))
dis.to_csv(f"{OUT}/44_PT_senmayo_by_disease.csv", encoding="utf-8-sig")
print("\nPT SenMayo by disease:"); print(dis.round(3).to_string())
pairs = [("Healthy_reference_tissue", "CKD"), ("Healthy_reference_tissue", "AKI"), ("Healthy_reference_tissue", "Deceased_donor")]
for a, b in pairs:
    u, p = stats.mannwhitneyu(pt.loc[pt.disease_cat == a, "senmayo"], pt.loc[pt.disease_cat == b, "senmayo"])
    print(f"  PT {a} vs {b}: P={p:.4g}")

# ---------- figures ----------
fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
top = ct.head(12).iloc[::-1]
axes[0].barh(top.index, top["senmayo_mean"], color="#4a7c96", alpha=0.85, xerr=top["senmayo_sd"] / np.sqrt(top["n"]))
axes[0].set_title("SenMayo senescence score by cell type\n(human kidney atlas, 304,989 nuclei)", fontsize=10.5)
axes[0].set_xlabel("Mean SenMayo score (z-based)"); axes[0].grid(alpha=0.3, axis="x")
axes[0].tick_params(axis="y", labelsize=8.5)
p16t = ct.sort_values("p16_detect", ascending=False).head(12).iloc[::-1]
axes[1].barh(p16t.index, p16t["p16_detect"] * 100, color="#b04a4a", alpha=0.85)
axes[1].set_title("CDKN2A (p16) detection rate by cell type (%)", fontsize=10.5)
axes[1].set_xlabel("p16 detection (%)"); axes[1].grid(alpha=0.3, axis="x"); axes[1].tick_params(axis="y", labelsize=8.5)
fig.tight_layout()
fig.savefig(f"{OUT}/EN_figureB18_senmayo_celltype.png", dpi=300, bbox_inches="tight")

fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))
ax = axes[0]
ax.scatter(donor_lvl["age"], donor_lvl["senmayo"], s=donor_lvl["n_cells"] / 8, c="#4a7c96", alpha=0.75, edgecolor="white", lw=0.5)
k = np.polyfit(donor_lvl["age"], donor_lvl["senmayo"], 1)
xs = np.linspace(donor_lvl["age"].min(), donor_lvl["age"].max(), 50)
ax.plot(xs, np.poly1d(k)(xs), "--", color="#555555", lw=1.4)
ax.set_title(f"Healthy donors: PT SenMayo vs donor age\n(n={len(donor_lvl)} donors, Spearman r={r_h:.3f}, P={p_h:.3g})", fontsize=10.5)
ax.set_xlabel("Donor age (years)"); ax.set_ylabel("PT-cell SenMayo score"); ax.grid(alpha=0.3)
ax = axes[1]
order = [d for d in ["Healthy_reference_tissue", "Healthy_stone_donor", "Tumor_nephrectomy", "Diabetes_mellitus_resilient", "CKD", "AKI", "Deceased_donor"] if d in set(pt.disease_cat)]
data = [pt.loc[pt.disease_cat == d, "senmayo"] for d in order]
bp = ax.boxplot(data, labels=[d.replace("_", "\n") for d in order], showfliers=False, patch_artist=True)
for patch, c in zip(bp["boxes"], ["#4a7c96"] * 2 + ["#7c8f4a"] * 3 + ["#b04a4a"] * 2):
    patch.set_facecolor(c); patch.set_alpha(0.6)
ax.set_title("PT-cell SenMayo score by disease category", fontsize=10.5)
ax.set_ylabel("PT-cell SenMayo score"); ax.tick_params(axis="x", labelsize=7.5); ax.grid(alpha=0.3, axis="y")
fig.tight_layout()
fig.savefig(f"{OUT}/EN_figureB19_PT_senmayo_age.png", dpi=300, bbox_inches="tight")
print("saved figures B18, B19 + CSVs")
