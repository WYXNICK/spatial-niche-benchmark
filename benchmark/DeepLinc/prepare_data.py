
import scanpy as sc
import pandas as pd
import numpy as np
import os
import argparse
import scipy.sparse as sp
from scipy.spatial.distance import pdist, squareform

def parse_args():
    parser = argparse.ArgumentParser(description="Prepare data for DeepLinc from H5AD")
    parser.add_argument('--input', '-i', type=str, required=True, help='Path to input .h5ad file')
    parser.add_argument('--outdir', '-o', type=str, required=True, help='Output directory for CSV files')
    parser.add_argument('--cell_type_key', '-k', type=str, default='cell_type', help='Column name for cell types in obs')
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Check paths
    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Input file not found: {args.input}")
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    print(f"Loading data from {args.input}...")
    adata = sc.read_h5ad(args.input)
    print(f"Data loaded: {adata.shape}")

    # --- 1. Basic Cleaning for Real Data ---
    # Convert gene names to unique if duplicates exist
    adata.var_names_make_unique()
    
    # 1. Prepare Counts (counts.csv)
    print("Generating counts.csv...")
    if sp.issparse(adata.X):
        # Determine chunk size to avoid memory error if matrix is huge
        # For 20k x 6k, dense is ~960MB (float64) or 480MB (float32), manageable.
        X = adata.X.toarray()
    else:
        X = adata.X
    
    # Save using float32 to save disk space and time if needed, but CSV is text.
    # Just standard save.
    counts_df = pd.DataFrame(X, columns=adata.var_names)
    counts_path = os.path.join(args.outdir, "counts.csv")
    counts_df.to_csv(counts_path, index=False)
    print(f"Saved {counts_path}")
    
    # Release memory
    del X, counts_df 

    # 2. Prepare Coordinates (coord.csv)
    print("Generating coord.csv...")
    if 'spatial' in adata.obsm:
        coords = adata.obsm['spatial']
    elif 'x' in adata.obs and 'y' in adata.obs:
        coords = adata.obs[['x', 'y']].values
    else:
        raise ValueError("Cannot find spatial coordinates in .obsm['spatial'] or .obs[['x', 'y']]")

    coord_df = pd.DataFrame(coords, columns=['x', 'y'])
    coord_path = os.path.join(args.outdir, "coord.csv")
    coord_df.to_csv(coord_path, index=False)
    print(f"Saved {coord_path}")

    # 3. Prepare Adjacency Matrix (adj.csv) - Memory Optimized
    print("Generating adj.csv (this may take a while for large datasets)...")
    
    # Optimization: Use float32 to save memory for the dense distance matrix
    # 20k * 20k * 4 bytes = ~1.6 GB. If float64, it's ~3.2 GB.
    dist_matrix = squareform(pdist(coords, 'euclidean')).astype(np.float32)

    link_num = 10
    dist_cutoff_arg = 0.95

    # Vectorized logic where possible, but loop is needed for KNN selection logic as per DeepLinc
    dist_neighbor_n = []
    
    # To speed up, we can use partition instead of full sort
    print("  Calculating distance cutoff...")
    for i in range(dist_matrix.shape[0]):
        # We need the distance to the 3rd nearest neighbor (excluding self)
        # Using partition is O(N) vs Sort O(NlogN)
        row = dist_matrix[i, :]
        # partition puts the k-th smallest element at position k
        # We want link_num-th element (0 is self, 1 is 1st NN, etc)
        # link_num=3 means we look at index 3 (4th element, 0-based: 0,1,2,3)
        # DeepLinc logic: neighbor_n_index list logic is complex, simplifying to exact equivalent:
        pass 
        # Actually, let's stick to sort for correctness unless it's too slow. 20k sorts is fine (~10s).
    
    # Optimized loop
    for i in range(dist_matrix.shape[0]):
        tmp = dist_matrix[i, :]
        # Select link_num-th nearest neighbor distance
        # sort [0.0, d1, d2, d3...], we want d(link_num)
        k_val = np.partition(tmp, link_num)[link_num]
        
        # Get neighbors
        neighbor_indices = np.where(tmp <= k_val)[0]
        # Remove self
        neighbor_indices = neighbor_indices[neighbor_indices != i]
        
        dist_neighbor_n.extend(dist_matrix[i, neighbor_indices].tolist())

    cutoff_distance = np.percentile(dist_neighbor_n, dist_cutoff_arg * 100)
    print(f"  Distance cutoff: {cutoff_distance:.2f}")

    # Create Adjacency Matrix
    print("  Building adjacency matrix...")
    # Condition 1: Distance < cutoff
    adj_1 = (dist_matrix < cutoff_distance).astype(np.int8) # 0 or 1
    np.fill_diagonal(adj_1, 0)
    
    # Condition 2: Is KNN
    # We iterate again to build adj_2 to save memory instead of making a full new float matrix
    # Or purely numpy:
    # Need to find k-th value for each row
    k_vals = np.partition(dist_matrix, link_num, axis=1)[:, link_num]
    # Broadcast comparison
    adj_2 = (dist_matrix <= k_vals[:, None]).astype(np.int8)
    np.fill_diagonal(adj_2, 0)

    # Combine (AND logic)
    adj_final = (adj_1 + adj_2) == 2
    adj_final = adj_final.astype(np.int8)
    
    # Symmetrize
    adj_final = ((adj_final + adj_final.T) > 0).astype(int)

    adj_df = pd.DataFrame(adj_final)
    adj_path = os.path.join(args.outdir, "adj.csv")
    # Write in chunks or directly if memory allows
    adj_df.to_csv(adj_path, index=False)
    print(f"Saved {adj_path}")
    
    del dist_matrix, adj_1, adj_2, adj_final

    # 4. Prepare Cell Labels (cell_type.csv)
    print("Generating cell_type.csv...")
    if args.cell_type_key in adata.obs:
        cell_types = adata.obs[args.cell_type_key].values
    else:
        print(f"Warning: '{args.cell_type_key}' column not found. Using 'Unknown'.")
        cell_types = ['Unknown'] * adata.shape[0]

    unique_types = np.unique(cell_types)
    type_to_id = {t: i for i, t in enumerate(unique_types)}
    class_ids = [type_to_id[t] for t in cell_types]

    cell_label_df = pd.DataFrame({
        'Cell_ID': adata.obs_names,
        'Cell_class_id': class_ids,
        'Cell_class_name': cell_types
    })
    cell_type_path = os.path.join(args.outdir, "cell_type.csv")
    cell_label_df.to_csv(cell_type_path, index=False)
    print(f"Saved {cell_type_path}")

    print("Data preparation complete.")

if __name__ == "__main__":
    main()
