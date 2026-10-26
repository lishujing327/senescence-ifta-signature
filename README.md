# A Senescence-Related Gene Signature for the Molecular Reclassification of Inflamed Tubular Atrophy in Kidney Allograft Biopsies

Analysis code for the manuscript submitted to *Clinical Transplantation*
(Original Report; manuscript "EN_v2"). Starting from a curated library of
1,043 ageing-related genes (GenAge + CellAge), the pipeline identifies
ageing-related genes associated with inflammation in areas of tubular atrophy
(i-IFTA) in kidney allograft protocol biopsies, constructs a 5-gene signature
by machine learning, localizes it at single-nucleus resolution, validates it
in independent cohorts, tests its prognostic association with graft loss,
benchmarks the senescence prior against whole-genome and random feature sets,
and compares it with established transplant molecular classifiers.

All data are from public repositories (GEO and CZ CELLxGENE Discover); no
human subjects data are included in this repository.

**v1.1**: file and directory names converted to ASCII for robust uploading;
analysis logic, locked numbers and figure content are unchanged from v1.0.

**v1.2**: the three series-matrix files in `geo_data/` (GSE22459, GSE53605,
GSE7392; ~25 MB total) and the 22 MB derived expression matrix
`results/fibrosis_ARG/39_GSE22459_log2expr.csv` were removed to keep every
file comfortably below the GitHub web-upload size limit. The matrices are
public GEO data re-downloadable in minutes (see `geo_data/EXTERNAL_INPUTS.md`);
the expression matrix is regenerated automatically by `EN_39.py` at the start
of the pipeline. All code, locked result tables and figure sources are
unchanged from v1.1.

## Repository layout

```
senescence-ifta-signature/
├── EN_39.py                          Differential expression (GSE22459) + volcano/heatmap figures
├── EN_39b_replot.py                  PPI network and enrichment figures (offline re-plot from cached STRING output)
├── EN_40.py                          Machine-learning comparison, 5-gene signature, ROC/nomogram/calibration/DCA
├── EN_41.py                          Single-nucleus localization + immune deconvolution (needs atlas .h5ad)
├── EN_42.py                          External validation (GSE53605, GSE7392)
├── 46_GSE21374_signature_outcome_EN.py   Prognostic association with graft loss (GSE21374)
├── 47_prior_value_benchmark_EN.py        Benchmarking of the senescence prior (whole-genome top-k vs random null)
├── 48_classifier_overlap_EN.py           Overlap with five established transplant molecular classifiers
├── daimon_runtime.py                 Font-configuration shim (no-op outside the Kimi Work runtime)
├── geo_data/                         Series matrices and probe mapping (small inputs; see EXTERNAL_INPUTS.md)
├── results/fibrosis_ARG/             Locked intermediate results (CSV) and English figure sources (PNG)
│   └── 49_figure_assembly_v2.py      Assembly of main/supplementary figures from figure sources
└── LICENSE                           MIT
```

## Environment

- Python ≥ 3.10
- Packages: `numpy`, `pandas`, `scipy`, `scikit-learn`, `matplotlib`, `h5py`
- The scripts were developed and run on Windows; paths are relative — **run every
  script from the repository root**.

```bash
pip install numpy pandas scipy scikit-learn matplotlib h5py
```

## Reproduction order

```bash
python EN_39.py                          # DEGs, ARDEGs, Figure 1 sources
python EN_39b_replot.py                  # Supplementary Figures S1-S2 sources (offline)
python EN_40.py                          # Models, signature, Figure 2 sources, Table S2
python EN_41.py                          # Figure 3 sources (requires atlas .h5ad; see EXTERNAL_INPUTS.md)
python EN_42.py                          # External validation, Supplementary Figure S5, Table S3
python 46_GSE21374_signature_outcome_EN.py  # Prognostic analysis, Main Figure 4 source
python 47_prior_value_benchmark_EN.py       # Prior benchmarking, Supplementary Figure S6
python 48_classifier_overlap_EN.py          # Classifier orthogonality, Supplementary Figure S7, Table S4
python results/fibrosis_ARG/49_figure_assembly_v2.py  # Assemble main/supplementary figures
```

Every manuscript number can be verified directly from the locked CSV files in
`results/fibrosis_ARG/` without re-running any download or analysis.

## Key locked results (for quick verification)

| Item | Value | File |
|---|---|---|
| ARDEGs (i-IFTA vs Norm) | 21 | `39_ARDEGs.csv` |
| 5-fold CV AUC, 5-gene signature | 0.907 (4-gene 0.918) | `40_EN_model_AUC.csv`, `42_EN_validation_AUC.csv` |
| GSE53605 IF/TA vs normal AUC | 1.000 (n=10) | `42_EN_validation_AUC.csv` |
| GSE7392 subclinical fibrosis AUC | 0.629 (P=0.236) | `42_EN_validation_AUC.csv` |
| GSE21374 graft loss, HR per 1 SD | 1.57 (95% CI 1.19-2.07; P=0.0013) | `46_GSE21374_HR_summary.csv` |
| Prior vs whole-genome top-k AUC | 0.896 vs 0.772-0.822 | `47_prior_benchmark_summary.csv` |
| Prior vs random 21-gene null | 100th percentile | `47_random_null_21genes.csv` |
| Overlap with 5 established classifiers | 0 genes (all panels) | `48_classifier_overlap_summary.csv` |

## Notes and caveats

- The prior-benchmarking analysis (script 47) is asymmetric by construction:
  the 21 ARDEGs were screened on all 41 modelling biopsies, whereas
  whole-genome features were selected within cross-validation folds. It
  quantifies the information content of the curated ageing library rather
  than an unbiased head-to-head comparison (see manuscript §2.10).
- For the Einecke 2010 classifier, the published 30 most-frequently-used
  genes (its Table 3) were used to represent the classifier gene list.
- Cox regression and log-rank tests are implemented directly with `scipy`;
  no external survival package is required.

## Citation

If you use this code, please cite the associated manuscript (in submission,
*Clinical Transplantation*) and this repository.

## License

MIT (see LICENSE).
