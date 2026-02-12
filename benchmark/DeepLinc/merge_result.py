
import scanpy as sc
import pandas as pd
import os

# Paths
base_dir = "../DeepLinc"
h5ad_path = os.path.join(base_dir, "dataset", "lymph_node_niche", "lymph_node_niche_annotated.h5ad")
result_dir = os.path.join(base_dir, "lymph_node_niche")
label_csv_path = os.path.join(result_dir, "label.csv")
output_h5ad_path = os.path.join(result_dir, "DeepLinc_lymph_node_niche.h5ad")

print(f"Loading data from {h5ad_path}...")
adata = sc.read_h5ad(h5ad_path)
print(f"Data loaded: {adata.shape}")

print(f"Loading labels from {label_csv_path}...")
if not os.path.exists(label_csv_path):
    raise FileNotFoundError(f"Label file not found: {label_csv_path}")

# Load the label CSV
# label.csv format:
# cell_id,cluster_id
# 0,0
# 1,1
# ...
label_df = pd.read_csv(label_csv_path)

# Ensure the cell_id matches the adata.obs index order
# Since the DeepLinc code (prepare_data.py) assumed adata.obs order was preserved, 
# and DeepLinc usually outputs in the same order as input (or with index), we need to be careful.
# DeepLinc's label.csv has 'cell_id' column which seems to be a 0-based integer index.
# We should map these 0-based indices back to adata.obs_names.

# Create a mapping from integer index to cluster_id
# Assuming label_csv 'cell_id' corresponds to the integer position in the original input
label_df['cell_id'] = label_df['cell_id'].astype(int)
label_df.set_index('cell_id', inplace=True)

# Sort by index just in case
label_df.sort_index(inplace=True)

# Assign to adata.obs
# We trust that the 0-th row in label.csv corresponds to the 0-th cell in adata.obs
# Verify lengths match
if len(label_df) != adata.n_obs:
    print(f"Warning: Label file has {len(label_df)} rows, but AnnData has {adata.n_obs} cells.")
    # If lengths differ, we might need to match by something else if possible, 
    # but DeepLinc usually works on the full dataset passed to it.

# Create 'pred_niche' column
# Initialize with a default value or NaN
adata.obs['pred_niche'] = pd.NA

# Map values. Since label_df index is integer 0..N-1, we can just assign the column if sorted.
# Using .loc with implicit integer index of obs if we reset index, or just direct array assignment if aligned.
if len(label_df) == adata.n_obs:
    adata.obs['pred_niche'] = label_df['cluster_id'].values.astype(str)
else:
    # Only assign matching indices if lengths differ (unlikely if pipeline is consistent)
    # adata.obs['pred_niche'] is aligned with adata.obs_names
    # label_df is index 0..M
    # We assume adata was processed sequentially.
    pass

# Also ensure it is categorical
adata.obs['pred_niche'] = adata.obs['pred_niche'].astype('category')

print("Saving new h5ad...")
adata.write_h5ad(output_h5ad_path)
print(f"Saved to {output_h5ad_path}")
