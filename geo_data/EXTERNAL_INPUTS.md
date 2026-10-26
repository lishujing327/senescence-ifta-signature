# External inputs (download required)

Four input files are **not** included in this repository because of their size
(3–80 MB each; binary matrices do not belong in version control). All of them
are freely available from public repositories. Download each file and place it
in this directory (`geo_data/`) under the exact file name listed below; every
analysis script reads the files from there.

The small inputs that **are** included in this directory are:

- `asia_map.tsv` — GPL96/97/571 probe-to-gene symbol mapping used by
  `EN_39.py`, `EN_42.py` and `46_GSE21374_signature_outcome_EN.py`
- `ARG_list.txt` — curated ageing-related gene library (1,043 genes; GenAge + CellAge)
- `GenAge_symbols.txt`, `cellage3.tsv` — raw GenAge / CellAge source tables

## 1. GSE22459 series matrix (GSE22459_matrix.txt.gz, ~13 MB)

- Source: GEO accession [GSE22459](https://www.ncbi.nlm.nih.org/geo/query/acc.cgi?acc=GSE22459)
  (GPL571 platform; kidney allograft protocol biopsies, i-IFTA discovery cohort)
- Download: on the GEO page, scroll to **Download family** and fetch
  `GSE22459_series_matrix.txt.gz`, then rename it to
  `GSE22459_matrix.txt.gz` and place it in `geo_data/`
- Used by: `EN_39.py` (differential expression, ARDEG screening)

## 2. GSE53605 series matrix (GSE53605_matrix.txt.gz, ~6.5 MB)

- Source: GEO accession [GSE53605](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53605)
  (GPL10558 platform; IF/TA vs normal kidney)
- Download: `GSE53605_series_matrix.txt.gz` from **Download family**, rename to
  `GSE53605_matrix.txt.gz`, place in `geo_data/`
- Used by: `EN_42.py` (external validation cohort 1)

## 3. GSE7392 series matrix (GSE7392_matrix.txt.gz, ~6 MB)

- Source: GEO accession [GSE7392](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7392)
  (GPL96/97 platform; subclinical fibrosis in protocol biopsies)
- Download: `GSE7392_series_matrix.txt.gz` from **Download family**, rename to
  `GSE7392_matrix.txt.gz`, place in `geo_data/`
- Used by: `EN_42.py` (external validation cohort 2)

## 4. GSE21374 series matrix (GSE21374_series_matrix.txt.gz, ~80 MB)

- Source: GEO accession [GSE21374](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21374)
  (GPL96/97 platform; kidney transplant biopsies with graft outcome follow-up)
- Download: `GSE21374_series_matrix.txt.gz` from **Download family**, keep the
  file name `GSE21374_series_matrix.txt.gz`, place in `geo_data/`
- Used by: `46_GSE21374_signature_outcome_EN.py` (prognostic association with
  graft loss) and `47_prior_value_benchmark_EN.py`

## 5. Human kidney single-nucleus atlas (02d5db96-d67e-45f2-b556-b6bd99753665.h5ad, ~1 GB)

- Source: CZ CELLxGENE Discover — the *Human kidney* single-nucleus RNA-seq
  reference atlas (search "Human kidney" on
  [CELLxGENE Discover](https://cellxgene.cziscience.com/))
- Download: use **Download → .h5ad** on the dataset page, then rename the
  downloaded file to `02d5db96-d67e-45f2-b556-b6bd99753665.h5ad` (the local
  download identifier, as referenced in `EN_41.py`) and place it in
  `geo_data/`. Use the same atlas version as in the manuscript so that cell
  type labels match.
- Used by: `EN_41.py` only (single-nucleus localization of the 21 ARDEGs and
  immune deconvolution of the signature). All other scripts run without it.

## Note

Every manuscript number can be verified from the locked CSV files in
`results/fibrosis_ARG/` **without** downloading or re-running anything.
The downloads above are only needed to reproduce the analysis from raw data.
