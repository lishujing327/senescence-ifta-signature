"""Minimal shim for the Kimi Work `daimon_runtime` helper.

The analysis scripts call `setup_plot()` to configure matplotlib fonts.
In the original runtime it registers bundled CJK fonts; all figure labels in
this pipeline are in English, so a plain DejaVu configuration is sufficient.
Place this file at the repository root and run scripts from the root.
"""
import matplotlib


def setup_plot():
    matplotlib.rcParams.update({
        "figure.dpi": 110,
        "savefig.dpi": 150,
        "font.family": "DejaVu Sans",
        "axes.unicode_minus": False,
        "axes.linewidth": 0.8,
    })
