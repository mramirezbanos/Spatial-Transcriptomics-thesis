import os
import random
import time
import torch
import tracemalloc
import scanpy as sc
import pandas as pd
import seaborn as sns
import numpy as np
import argparse
from sklearn.decomposition import PCA
from pathlib import Path
import SEDR

os.environ['OMP_NUM_THREADS'] = '1'
os.environ['KMP_INIT_AT_FORK'] = 'FALSE'
os.environ['MKL_NUM_THREADS'] = '1'

# == Define main ==
def main(in_file, annotation_key, spatial_x, spatial_y, k_geom=15, resolutions=[0.4], pca_dims=[20], n_hvg=3000):
    
    # Track time and memory
    random_seed = 0
    SEDR.fix_seed = random_seed
    tracemalloc.start()
    t1 = time.time()

    # Paths
    in_path = Path("/faststorage/project/limitedAccessBladder/NOBACKUP/BIRC_thesis_students_spring26/marinarb/outputs/normalized_hvgs") / f"{in_file}.h5ad"
    out_path = Path("/faststorage/project/limitedAccessBladder/NOBACKUP/BIRC_thesis_students_spring26/marinarb/outputs/SEDR/3000_hvg_Norm_50PC_celltype_r5_predclusters/")

    # == Load ==
    adata = sc.read_h5ad(in_path)
    print("\nLoaded:", in_path)
    adata.obs["ground_truth"] = adata.obs[annotation_key]
    adata = adata[~pd.isnull(adata.obs["ground_truth"])]
    n_cluster = len(set(adata.obs["ground_truth"])) # remove unassigned cells
    print(adata.shape)
    print(f"Number of ground truth clusters: {n_cluster}")
    print(f"Using {n_hvg} HVGs")

    # Update spatial coordinates
    if 'spatial' in adata.obs:
        del adata.obs['spatial']
    if 'spatial' in adata.obsm:
        del adata.obsm['spatial']
    adata.obsm['spatial'] = adata.obs[[spatial_x, spatial_y]].to_numpy()
    
    # == Required preprocessing ==
    # the data imported is already normalized by depth, log1p transformed and feature selected
    print("\nRunning gene-wise data scaling & PCA")
    adata.layers['count'] = adata.X.copy()
    sc.pp.scale(adata)
    
    adata_X = PCA(n_components=50, random_state=random_seed).fit_transform(adata.X) # sklearn PCA is used because PCA in scanpy is not stable.
    adata.obsm['X_pca'] = adata_X
    
    print("\nSolved. Continue...")

    # == Neighborhood graph ==
    print("\nInitialize SEDR spatial graph...")
    graph_dict = SEDR.graph_construction(
        adata,
        n=6)

    # == SEDR Training ==
    sedr_net = SEDR.Sedr(adata.obsm['X_pca'], graph_dict, mode='clustering', device='cpu')
    using_dec = True
    if using_dec:
        sedr_net.train_with_dec()
    else:
        sedr_net.train_without_dec()
    sedr_feat, _, _, _ = sedr_net.process()
    adata.obsm['SEDR'] = sedr_feat

    # == Clustering ==
    #os.environ['R_HOME'] = '/home/marinarb/miniforge3/envs/sedr_fix/lib/R'
    #SEDR.mclust_R(adata, n_cluster, use_rep='SEDR', key_added='SEDR')
    
    from sklearn.cluster import KMeans
    kmeans = KMeans(n_clusters=n_cluster, random_state=0, n_init=10)
    kmeans.fit(adata.obsm['SEDR'])
    adata.obs['SEDR'] = [f"SEDR_{int(x)}" for x in kmeans.labels_]
    
    # == Save results ==
    adata.obs['predict_cluster'] = adata.obs['SEDR'].astype('category')

    obs_csv_path = out_path / f"{in_file}_SEDR_predclusters.csv"

    adata.obs[['predict_cluster']].to_csv(obs_csv_path, index=True)

    print(f"Per-cell cluster CSV saved: {obs_csv_path}")
    print(f"\nResults for {in_file} saved.")

    # Measure peak memory usage and execution time
    t2 = time.time()
    t = t2 - t1
    size, peak = tracemalloc.get_traced_memory()
    memory = peak / 1024 / 1024
    tracemalloc.stop()
    print("\nExecution time:", t)
    print("\nPeak memory usage:", memory)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_file', required=True)
    parser.add_argument('--annotation_key', required=True)
    parser.add_argument('--spatial_x', required=True)
    parser.add_argument('--spatial_y', required=True)
    parser.add_argument('--k_geom', type=int, default=15)
    parser.add_argument('--resolutions', nargs='+', type=float, default=[0.76])
    parser.add_argument('--pca_dims', nargs='+', type=int, default=[20])
    parser.add_argument('--n_hvg', type=int, default=3000)

    args = parser.parse_args()
    main(args.in_file, args.annotation_key, args.spatial_x, args.spatial_y, args.k_geom, args.resolutions, args.pca_dims, args.n_hvg)
