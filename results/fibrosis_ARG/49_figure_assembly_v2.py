# -*- coding: utf-8 -*-
"""English main/supplementary figure assembly v2 (A1-A3 + A6/STARD integrated):
Main Figure 4 = GSE21374 outcome association (former external validation demoted to FigS5; FigS6 prior benchmark, FigS7 classifier orthogonality added).
Main Figures 1-3 identical to v1."""
import os, shutil
from PIL import Image, ImageDraw, ImageFont
import matplotlib

FONT = ImageFont.truetype(
    os.path.join(os.path.dirname(matplotlib.__file__), "mpl-data", "fonts", "ttf", "DejaVuSans-Bold.ttf"), 48)

SRC = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(SRC, "figures_english")
SUP = os.path.join(OUT, "Supplementary")
os.makedirs(SUP, exist_ok=True)

LBL, MARG, TOP = 48, 60, 95

def load(name): return Image.open(os.path.join(SRC, name)).convert("RGB")
def scale_h(im, h): return im.resize((round(im.width * h / im.height), h), Image.LANCZOS)
def scale_w(im, w): return im.resize((w, round(im.height * w / im.width)), Image.LANCZOS)

def panel(canvas, im, x, y, label):
    canvas.paste(im, (x, y))
    ImageDraw.Draw(canvas).text((x, y - LBL - 14), label, fill="black", font=FONT)

def new_canvas(w, h): return Image.new("RGB", (w, h), "white")

# Main Figure 1: a volcano + b heatmap (identical to v1)
b1, b2 = scale_h(load("EN_figureB1_volcano.png"), 1100), scale_h(load("EN_figureB2_heatmap.png"), 1100)
W = MARG + b1.width + 100 + b2.width + MARG
c = new_canvas(W, TOP + 1100 + 40)
panel(c, b1, MARG, TOP, "a"); panel(c, b2, MARG + b1.width + 100, TOP, "b")
c.save(os.path.join(OUT, "Main_Figure1.png"), dpi=(300, 300)); print("MF1", c.size)

# Main Figure 2: a models + b ROC + c distribution (identical to v1)
bot_b = scale_h(load("EN_figureB6_ROC.png"), 1000)
bot_c = scale_h(load("EN_figureB7_signature_distribution.png"), 1000)
row2_w = bot_b.width + 100 + bot_c.width
top_im = scale_w(load("EN_figureB5_model_comparison.png"), row2_w)
W = MARG + row2_w + MARG
y2 = TOP + top_im.height + 90
c = new_canvas(W, y2 + 1000 + 40)
panel(c, top_im, MARG, TOP, "a")
panel(c, bot_b, MARG, y2, "b"); panel(c, bot_c, MARG + bot_b.width + 100, y2, "c")
c.save(os.path.join(OUT, "Main_Figure2.png"), dpi=(300, 300)); print("MF2", c.size)

# Main Figure 3: a cell localization + b immune deconvolution (identical to v1)
a = scale_h(load("EN_figureB10_cell_localization.png"), 1100)
b = scale_h(load("EN_figureB11_immune_deconvolution.png"), 1100)
W = MARG + a.width + 80 + b.width + MARG
c = new_canvas(W, TOP + 1100 + 40)
panel(c, a, MARG, TOP, "a"); panel(c, b, MARG + a.width + 80, TOP, "b")
c.save(os.path.join(OUT, "Main_Figure3.png"), dpi=(300, 300)); print("MF3", c.size)

# Main Figure 4 (added in v2): GSE21374 outcome association (source figure carries a/b labels)
shutil.copy(os.path.join(SRC, "EN_figureB13_outcome_GSE21374.png"), os.path.join(OUT, "Main_Figure4.png"))
print("MF4 = outcome (B13)")

# Supplementary: S1-S4 unchanged; S5 external validation; S6 prior benchmark; S7 classifier orthogonality; S8 donor age (A6); S9 STARD flow diagram; S10 SenMayo by cell type (A7); S11 PT senescence age/disease (A7)
for src, dst in [
    ("EN_figureB3_PPI_network.png", "FigS1_PPI_network.png"),
    ("EN_figureB4_enrichment.png", "FigS2_Functional_enrichment.png"),
    ("EN_figureB8_nomogram.png", "FigS3_Nomogram.png"),
    ("EN_figureB9_calibration_DCA.png", "FigS4_Calibration_DCA.png"),
    ("EN_figureB12_external_validation.png", "FigS5_External_validation.png"),
    ("EN_figureB14_prior_value.png", "FigS6_Prior_value_benchmark.png"),
    ("EN_figureB15_classifier_overlap.png", "FigS7_Classifier_orthogonality.png"),
    ("EN_figureB17_donor_age_GSE147451.png", "FigS8_Donor_age_correlation.png"),
    ("EN_figureB16_STARD_flow.png", "FigS9_STARD_flow_diagram.png"),
    ("EN_figureB18_senmayo_celltype.png", "FigS10_SenMayo_by_celltype.png"),
    ("EN_figureB19_PT_senmayo_age.png", "FigS11_PT_senescence_age_disease.png"),
]:
    shutil.copy(os.path.join(SRC, src), os.path.join(SUP, dst))
print("Supplementary 11 done")

# Individual figure sources
IND = os.path.join(OUT, "Individual_figures")
os.makedirs(IND, exist_ok=True)
for f in sorted(os.listdir(SRC)):
    if f.startswith("EN_figureB") and f.endswith(".png"):
        shutil.copy(os.path.join(SRC, f), os.path.join(IND, f))
print("Individual EN figures backed up")
print("OUT:", OUT)
