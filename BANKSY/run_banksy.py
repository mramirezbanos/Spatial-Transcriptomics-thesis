import os, csv
from pathlib import Path
import gc
import sys
import time
import tracemalloc
import random

import numpy as np
import pandas as pd

import scipy.sparse as sp
import scanpy as sc
import anndata as ad

import banksy.main as bm
from banksy.initialize_banksy import initialize_banksy
from banksy.run_banksy import run_banksy_multiparam
#from banksy_utils.color_lists import spagcn_color
import matplotlib.pyplot as plt

from sklearn.metrics.cluster import contingency_matrix
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, homogeneity_completeness_v_measure, fowlkes_mallows_score
from sklearn.preprocessing import StandardScaler

# == Define helpers ==
def params_to_string(lambda_list, k_geom, resolutions, pca_dims, n_hvg):
    lam_str = '-'.join([f"{lam:.1f}".replace('.', 'p') for lam in lambda_list])
    res_str = '-'.join([f"{res:.1f}".replace('.', 'p') for res in resolutions])
    pca_str = '-'.join(map(str, pca_dims))
    return f"hvg{n_hvg}_{lam_str}_k{k_geom}_res{res_str}_pca{pca_str}"

def calculate_purity(ground_truth, cluster):
    cm = contingency_matrix(ground_truth, cluster)
    return np.sum(np.amax(cm, axis=0)) / np.sum(cm)

def calculate_v_measure(ground_truth, cluster):
    homogeneity, completeness, vmeasure = homogeneity_completeness_v_measure(ground_truth, cluster)
    return homogeneity, completeness, vmeasure
    
# == Define main ==
def main(in_file, lambda_list=[0.2, 0.5], k_geom=15, resolutions=[0.4, 0.6], pca_dims=[20]):

    # Track time and memory
    tracemalloc.start()
    t1 = time.time()

    # Paths
    in_path = Path("/faststorage/project/limitedAccessBladder/NOBACKUP/BIRC_thesis_students_spring26/marinarb/outputs/normalized_hvgs") / f"{in_file}.h5ad"
    out_path = Path("/faststorage/project/limitedAccessBladder/NOBACKUP/BIRC_thesis_students_spring26/marinarb/outputs/BANKSY/not_log1p_3000_hvg_gt_celltype_r5/predicted_clusters/")

    # == Load ==
    adata = sc.read_h5ad(in_path)
    print("\nLoaded:", in_path)
    print(adata.shape)
    
    # Sanity: ensure log1p transformed data and feature selection done
    print("\nRunning sanity checks on data format...")
    #if "log1p" not in adata.uns or not np.issubdtype(adata.X.dtype, np.floating):
        #raise ValueError("Missing log-transfromation and/or X not float)")
    
    if "highly_variable" not in adata.var.columns:
        raise ValueError("Missing adata.var['highly_variable']")
    
    if "highly_variable" in adata.var.columns:
        hvg_mask = adata.var["highly_variable"].to_numpy()
        n_hvg = int(hvg_mask.sum())
        if n_hvg == 0:
            raise ValueError("No valid HVGs found")
    
    # Sanity: check for required fields
    required = ["CenterX_global_px", "CenterY_global_px", "core_key"]
    
    if not all(col in adata.obs.columns for col in required):
        missing = [col for col in required if col not in adata.obs.columns]
        raise KeyError(f"Missing required obs columns: {missing}")
    
    print("\nSanity checks done")
    gc.collect()

    # == BANKSY hyperparameters ==
    # fixed
    nbr_weight_decay = 'scaled_gaussian'
    core_spacing = 1e7
    cluster_algorithm = 'leiden'
    n_clust = 0
    
    # clusters (equal to ground truth)
    annotation_key = 'ground_truth'
    n_cluster = len(set(adata.obs["ground_truth"]))
    
    print(f"Configuration used: lambda {lambda_list} | neighbours {k_geom}, | leiden clustering resolution {resolutions} | PCs {pca_dims}. Using {n_hvg} HVGs")

    # == Prepare BANKSY on HVGs (working copy) =
    print("\nPreparing data subset of HVGs for BANKSY...")
    adata.layers['counts'] = adata.X.copy()
    adata_b = adata[:, hvg_mask].copy()
    print("\nBANKSY using (HVG subset):", adata_b.shape)
    
    # -- Z-score matrix transformation (gene-wise) --
    print("\nTranform matrix using gene-wise z-score (dense matrix; it can take long time)...")
    Z = bm.zscore(adata_b.X, axis=0).astype(np.float32, copy=False)
    adata_b.X = Z # work with z-scores
    del Z
    print("Z-scored X dtype:", adata_b.X.dtype, "; Shape:", adata_b.X.shape)
    print("\nMatrix succesfully transformed!")
    
    gc.collect()
    
    # -- Prevent cross-core neighbours --
    # Sanity: safe numbers
    print("\nPeventing cross-core neighbours...")
    coords_x = pd.to_numeric(adata_b.obs["CenterX_global_px"], errors='coerce').to_numpy()
    coords_y = pd.to_numeric(adata_b.obs["CenterY_global_px"], errors='coerce').to_numpy()
    
    if np.isnan(coords_x).any() or np.isnan(coords_y).any():
        raise ValueError("Invalid coordinates (NaN detected)")
    
    # Adding spatial separation between cores
    if "core_codes" not in adata_b.obs.columns:
        
        print("\nCores not labeled. Generating core_ids...")
        core_ids = pd.Categorical(adata_b.obs["core_key"]).codes
        
        print("\nAdding spacing...")
        spatial_x = coords_x.astype(np.float64)  + core_ids.astype(np.float64) * core_spacing
        spatial_y = coords_y.astype(np.float64)
    
        # Save modifications
        adata_b.obs["x_banksy"] = spatial_x
        adata_b.obs["y_banksy"] = spatial_y
        adata_b.obsm["spatial_banksy"] = np.column_stack([spatial_x, spatial_y])
        coord_keys = ("x_banksy", "y_banksy", "spatial_banksy")
        print("\nCoordinate keys updated to:", coord_keys, "in HVG susbet")
        
    else: 
        if "x_banksy" and "y_banksy" in adata_b.obs.columns:
            print("\nCores were aleady separated")
    print("\nSolved. Continue...")
   
    ## == BANSKY spatial graph ==
    print("\nInitialize BANKSY spatial graph...")
    banksy_dict = initialize_banksy(
        adata_b,
        coord_keys=coord_keys,
        num_neighbours=k_geom,
        nbr_weight_decay=nbr_weight_decay,
        plt_edge_hist=False,
        plt_nbr_weights=False,
        plt_agf_angles=False,
        plt_theta=False
    )
    print("\nDone. Continue...")
    
    gc.collect()

    # == BANKSY clustering ==
    print("\nStarting BANKSY clustering and dimensionality reduction. Using Leiden as cluster algorithm")
    
    colorlist = plt.cm.tab20(np.linspace(0, 1, 100)).tolist()
    
    results_df = run_banksy_multiparam(
        adata_b,
        banksy_dict,
        lambda_list,
        resolutions,
        color_list = colorlist,
        filepath = out_path,
        max_m = 1,
        key = coord_keys,
        pca_dims = pca_dims,
        annotation_key = annotation_key,
        max_labels = 100,
        cluster_algorithm = cluster_algorithm,
        match_labels = False,
        
        savefig = False,
        add_nonspatial = False,
        variance_balance = False,
    )
    print("\nBANKSY completed. Proceed to save results as .csv")
    gc.collect()
       
    
    # Save results summary
    results_summary = []
    for config_name in results_df.index:
        
        print(f"Processing {config_name}...")
        
        # Extract configuration used
        adata_config = results_df.loc[config_name, 'adata']
        adata_config.obs['Banksy'] = adata_config.obs[config_name].astype('category')
        
        # Save predicted cluster for each cell
        adata_config.obs['predict_cluster'] = adata_config.obs['Banksy'].astype('category')
        
        obs_csv_path = out_path / f"{in_file}_{config_name}_predclusters.csv"
        
        adata_config.obs[['predict_cluster']].to_csv(obs_csv_path, index=True)
        
        print(f"Per-cell cluster CSV saved: {obs_csv_path}")
            
        # Obtain performance metrics
        gt = adata_config.obs['ground_truth']
        pred = adata_config.obs['Banksy']
    
        n_clusters_pred = len(np.unique(pred))
        n_clusters_gt = len(np.unique(gt))    
    
        print(f"Obtaining accuracy metrics for {config_name}...")
        ari = adjusted_rand_score(gt, pred)
        nmi = normalized_mutual_info_score(gt, pred)
        fmi = fowlkes_mallows_score(gt, pred)
        purity = calculate_purity(gt, pred)
        homog, comp, vmeas = calculate_v_measure(gt, pred)
                
        # Save .csv
        print(f"Saving results for {config_name} ...")
        results_summary.append({
            'file_name': f"{in_file}_{config_name}_banksy_summary.csv", 
            'config': config_name,
            'ARI': round(ari, 4),
            'NMI': round(nmi, 4),
            'FMI': round(fmi, 4),
            'Purity': round(purity, 4),
            'Homogeneity': round(homog, 4),
            'Completeness': round(comp, 4),
            'Vmeasure': round(vmeas, 4),
            'n_clusters_pred': n_clusters_pred,
            'n_clusters_gt': n_clusters_gt,
            'n_cells': len(adata_config)
        })
        
        print(f"ARI {config_name}: {ari:.4f}")
        
    
    param_str = params_to_string(lambda_list, k_geom, resolutions, pca_dims, n_hvg)
    out_file = f"{in_file}_{param_str}.csv"
    
    csv_path = out_path / out_file
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    
    summary_df = pd.DataFrame(results_summary)
    summary_df.to_csv(csv_path, index=False)
    
    print(f"\nResults for {in_file} saved in:")
    print(f" Summary: {csv_path}")
    print(f"  Per config: {out_file}")
    
    
    # Measure peak memory usage and execution time
    t2 = time.time()
    t = t2 - t1
    size, peak = tracemalloc.get_traced_memory()
    memory = peak / 1024 / 1024
    tracemalloc.stop()
    print("\nExecution time:", t)
    print("\nPeak memory usage:", memory)
    
import argparse
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_file', required=True)
    parser.add_argument('--lambda_list', nargs='+', type=float, default=[0.2, 0.5])
    parser.add_argument('--k_geom', type=int, default=15)
    parser.add_argument('--resolutions', nargs='+', type=float, default=[0.4, 0.6])
    parser.add_argument('--pca_dims', nargs='+', type=int, default=[20])
    args = parser.parse_args()
    main(args.in_file, args.lambda_list, args.k_geom, args.resolutions, args.pca_dims)

