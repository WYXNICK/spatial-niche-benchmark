# Annotations

This folder contains cell-level annotations used in the benchmark.

## Files
- `lymph_node_annotations.tsv` (recommended; tab-separated)
- `lymph_node_annotations.csv` (CSV copy)

## Table schema
Each row corresponds to one cell.

Columns:
- `cell_id`: cell identifier (matches `adata.obs_names` in the benchmark-ready `.h5ad`)
- `x`, `y`: spatial coordinates from `adata.obs["x"]` and `adata.obs["y"]`
- `cell_type_annotation`: cell-type label
- `niche_annotation`: niche label (4 classes: B cell Zone, T cell Zone, Medulla, Germinal Center)

## Load and merge (example)
```python
import pandas as pd

ann = pd.read_csv("annotations/lymph_node_annotations.tsv", sep="\t").set_index("cell_id")

# adata.obs_names should match ann.index (cast to str if needed)
adata.obs["cell_type_annotation"] = ann.loc[adata.obs_names.astype(str), "cell_type_annotation"].values
adata.obs["niche_annotation"] = ann.loc[adata.obs_names.astype(str), "niche_annotation"].values