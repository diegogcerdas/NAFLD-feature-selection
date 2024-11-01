# Feature Selection on Gene Expression Data - NAFLD Biomarker Identification
Repo owner: Diego García Cerdas.

This repo performs feature selection on gene expression datasets to identify biomarkers for Nonalcoholic Fatty Liver Disease (NAFLD). Specifically, it looks at the following datasets: *GSE49541*, *GSE48452*, *GSE89632*, and *GSE83452*.

## Repo Structure

### Data Exploration

Jupyter notebooks in `data_exploration/` go through each dataset to extract relevant information: (1) gene expression intensity values, (2) the condition of each patient, and (3) a symbol for each gene in the dataset.
- The resulting preprocessed datasets can be loaded through methods in `datasets.py` (e.g., `load_GSE49541()`).
- The file `utils.py` contains a method for converting GenBank accession IDs to gene symbols.

### Feature selection

The notebook `feature_selection.ipynb` identifies relevant genes for NAFL (steatosis) vs NASH, as well as Control vs NASH, using Recursive Feature Elimination (RFE). It also identifies biological pathways that are relevant across datasets.