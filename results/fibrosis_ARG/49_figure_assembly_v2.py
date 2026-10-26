# -*- coding: utf-8 -*-
"""英文版投稿主图拼版 v2（A1-A3 整合版）：
Main Figure 4 = GSE21374 结局关联（原外部验证降为补充 FigS5；新增 FigS6 先验对照、FigS7 分类器正交）。
Main Figure 1-3 与 v1 完全一致。"""
import os, shutil
from PIL import Image, ImageDraw, ImageFont
import matplotlib

FONT = ImageFont.truetype(
    os.path.join(os.path.dirname(matplotlib.__file__), "mpl-data", "fonts", "ttf", "DejaVuSans-Bold.ttf"), 48)

SRC = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(SRC, "投稿图片_英文版_20261024")
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

# Main Figure 1: a volcano + b heatmap（与 v1 一致）
b1, b2 = scale_h(load("EN_figureB1_volcano.png"), 1100), scale_h(load("EN_figureB2_heatmap.png"), 1100)
W = MARG + b1.width + 100 + b2.width + MARG
c = new_canvas(W, TOP + 1100 + 40)
panel(c, b1, MARG, TOP, "a"); panel(c, b2, MARG + b1.width + 100, TOP, "b")
c.save(os.path.join(OUT, "Main_Figure1.png"), dpi=(300, 300)); print("MF1", c.size)

# Main Figure 2: a models + b ROC + c distribution（与 v1 一致）
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

# Main Figure 3: a cell localization + b immune deconvolution（与 v1 一致）
a = scale_h(load("EN_figureB10_cell_localization.png"), 1100)
b = scale_h(load("EN_figureB11_immune_deconvolution.png"), 1100)
W = MARG + a.width + 80 + b.width + MARG
c = new_canvas(W, TOP + 1100 + 40)
panel(c, a, MARG, TOP, "a"); panel(c, b, MARG + a.width + 80, TOP, "b")
c.save(os.path.join(OUT, "Main_Figure3.png"), dpi=(300, 300)); print("MF3", c.size)

# Main Figure 4 (v2 新增): GSE21374 结局关联（图源自带 a/b 标注）
shutil.copy(os.path.join(SRC, "EN_figureB13_outcome_GSE21374.png"), os.path.join(OUT, "Main_Figure4.png"))
print("MF4 = outcome (B13)")

# Supplementary: S1-S4 不变；S5=外部验证（原主图4降级）；S6=先验对照；S7=分类器正交
for src, dst in [
    ("EN_figureB3_PPI_network.png", "FigS1_PPI_network.png"),
    ("EN_figureB4_enrichment.png", "FigS2_Functional_enrichment.png"),
    ("EN_figureB8_nomogram.png", "FigS3_Nomogram.png"),
    ("EN_figureB9_calibration_DCA.png", "FigS4_Calibration_DCA.png"),
    ("EN_figureB12_external_validation.png", "FigS5_External_validation.png"),
    ("EN_figureB14_prior_value.png", "FigS6_Prior_value_benchmark.png"),
    ("EN_figureB15_classifier_overlap.png", "FigS7_Classifier_orthogonality.png"),
]:
    shutil.copy(os.path.join(SRC, src), os.path.join(SUP, dst))
print("Supplementary 7 done")

# 单独图源
IND = os.path.join(OUT, "Individual_figures")
os.makedirs(IND, exist_ok=True)
for f in sorted(os.listdir(SRC)):
    if f.startswith("EN_figureB") and f.endswith(".png"):
        shutil.copy(os.path.join(SRC, f), os.path.join(IND, f))
print("Individual EN figures backed up")
print("OUT:", OUT)
