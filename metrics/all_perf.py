import os
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import anndata as ad
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import scipy.sparse as sp
from scipy.spatial import distance_matrix
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, homogeneity_completeness_v_measure, fowlkes_mallows_score, silhouette_score
from sklearn.metrics.cluster import contingency_matrix
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors

# == Define helpers ==
# Accuracy metrics 
def get_ari(gt, pred):
    ari = adjusted_rand_score(gt, pred)
    return ari

def get_nmi(gt, pred):
    nmi = normalized_mutual_info_score(gt, pred)
    return nmi
    
def get_fmi(gt, pred):
    fmi = fowlkes_mallows_score(gt, pred)
    return fmi

def get_purity(gt, pred):
    cm = contingency_matrix(gt, pred)
    return np.sum(np.amax(cm, axis=0)) / np.sum(cm)

def get_v_measure(gt, pred):
    homogeneity, completeness, vmeasure = homogeneity_completeness_v_measure(gt, pred)
    return homogeneity, completeness, vmeasure


# Spatial continuity metrics
def build_spatial_edges(coords, n_neighbors=10):
    coords = np.asarray(coords, dtype=float)
    n = coords.shape[0]
    if n < 2:
        return []
    k = min(n_neighbors + 1, n)
    nn = NearestNeighbors(n_neighbors=k)
    nn.fit(coords)
    _, indices = nn.kneighbors(coords)
    edges = set()
    for i in range(n):
        for j in indices[i, 1:]:
            a, b = sorted((i, j))
            edges.add((a, b))
    return list(edges)

def spatial_entropy_from_edges(labels, edges):
    labels = np.asarray(labels)
    if len(edges) == 0:
        return np.nan
    pair_counts = {}
    for i, j in edges:
        a = labels[i]
        b = labels[j]
        key = tuple(sorted((str(a), str(b))))
        pair_counts[key] = pair_counts.get(key, 0) + 1
    counts = np.array(list(pair_counts.values()), dtype=float)
    probs = counts / counts.sum()
    return -np.sum(probs * np.log(probs + 1e-12))
    
def nearest_neighbor_distance(point_idx, coords):
    pairwise_dist = distance_matrix(coords[[point_idx]], coords).ravel()
    pairwise_dist[point_idx] = np.inf
    return pairwise_dist.min()

def compute_spatial_pas(labels, coords, k_neighbors=10):
    labels = np.asarray(labels)
    coords = np.asarray(coords)
    
    misclassified_count = 0
    for i in range(len(labels)):
        pairwise_dist = distance_matrix(coords[[i]], coords).ravel()
        pairwise_dist[i] = np.inf
        neighbor_indices = np.argpartition(pairwise_dist, k_neighbors)[:k_neighbors]
        
        if np.sum(labels[neighbor_indices] != labels[i]) > (k_neighbors / 2):
            misclassified_count += 1
    
    return misclassified_count / len(labels)

def compute_spatial_chaos(labels, coords):
    labels = np.asarray(labels)
    coords = np.asarray(coords)
    
    scaler = StandardScaler()
    scaled_coords = scaler.fit_transform(coords)
    
    total_nn_distance = 0
    valid_clusters = 0
    
    for cluster_id in np.unique(labels):
        cluster_points = scaled_coords[labels == cluster_id]
        if len(cluster_points) > 2:
            cluster_nn_distances = [nearest_neighbor_distance(i, cluster_points) 
                                  for i in range(len(cluster_points))]
            total_nn_distance += np.sum(cluster_nn_distances)
            valid_clusters += 1
    
    return total_nn_distance / len(labels) if valid_clusters > 0 else 0

def compute_asw_spatial(adata, cluster_col, coords):
    dist_matrix = squareform(pdist(coords))
    cluster_labels = adata.obs[cluster_col].values
    
    return silhouette_score(dist_matrix, cluster_labels, metric='precomputed')

def compute_spatial_coherence_score(labels, coords, n_neighbors=10, n_permutations=100, random_state=0):
    labels = np.asarray(labels)
    coords = np.asarray(coords, dtype=float)

    if len(labels) != len(coords):
        raise ValueError("labels and coords must have the same length")

    if len(labels) < 3:
        return np.nan

    if len(np.unique(labels)) < 2:
        return np.nan

    edges = build_spatial_edges(coords, n_neighbors=n_neighbors)
    if len(edges) == 0:
        return np.nan

    observed_entropy = spatial_entropy_from_edges(labels, edges)

    rng = np.random.default_rng(random_state)
    permuted_entropies = []
    for _ in range(n_permutations):
        permuted_labels = rng.permutation(labels)
        ent = spatial_entropy_from_edges(permuted_labels, edges)
        permuted_entropies.append(ent)

    permuted_entropies = np.asarray(permuted_entropies, dtype=float)

    if np.std(permuted_entropies, ddof=1) == 0:
        return np.nan

    z_entropy = (observed_entropy - np.mean(permuted_entropies)) / np.std(permuted_entropies, ddof=1)
    scs = -z_entropy
    return scs

# == Define main function ==
results_summary = []

def main (in_file, out_path, meta_dir, file_predicted_clustering, method, annotation_key, spatial_x, spatial_y):
    # == Load data ==
    meta_file = f"{file_predicted_clustering}.csv"
    # in_path = Path("/faststorage/project/limitedAccessBladder/NOBACKUP/BIRC_thesis_students_spring26/marinarb/outputs/normalized_log1p_hvgs") / f"{in_file}.h5ad"
    out_path = Path(out_path)
    meta_dir = Path(meta_dir)
    meta_path = meta_dir / meta_file
    out_file = f"{file_predicted_clustering}_accuracy_summary.csv"
    csv_path = out_path / out_file
    
    # adata = sc.read_h5ad(in_path)
    # print("\nLoaded:", in_path)
    adata = in_file

   # -- Adding metadata to adata --
    df_meta = pd.read_csv(meta_path)
    df_meta = df_meta.set_index(df_meta.columns[0])
    
    adata = adata[adata.obs.index.isin(df_meta.index)].copy()
    df_meta = df_meta.reindex(adata.obs.index)
    
    # filter adata
    cells_to_keep = adata.obs.index.isin(df_meta.index)
    adata = adata[cells_to_keep].copy()
    adata.shape

    # sanity
    # if ground_truth in adata.obs:
    #     del adata.obs[ground_truth]
    
    # add predited clusters
    adata.obs["predict_cluster"] = df_meta["predict_cluster"].astype("category")
    adata.obs["ground_truth"] = adata.obs[annotation_key].astype("category")
    
    print("Total predicted clusters:", adata.obs['predict_cluster'].nunique())
    print("Total ground truth clusters:", adata.obs['ground_truth'].nunique())
    
    
    print("\nMetadata added to adata object:")
    print(adata.obs.head())

    # == Obtain accuracy and spatial continuity metrics ==    
    # Obtain performance metrics
    gt = adata.obs['ground_truth']
    pred = adata.obs['predict_cluster']
    
    n_clusters_pred = len(np.unique(pred))
    n_clusters_gt = len(np.unique(gt))    
    
    print(f"Obtaining accuracy metrics ...")
    ari = get_ari(gt, pred)
    nmi = get_nmi(gt, pred)
    fmi = get_fmi(gt, pred)
    purity = get_purity(gt, pred)
    homog, comp, vmeas = get_v_measure(gt, pred)
    
    print(ari, nmi, fmi, purity, homog, comp, vmeas)
    print(f"Done. Continue..")
        
    print(f"Computing spatial metrics ...")
    
    spatial_coords = np.column_stack([
        adata.obs[spatial_x].astype(float),
        adata.obs[spatial_y].astype(float)
    ])
    
    pas_score = compute_spatial_pas(pred.values, spatial_coords)
    chaos_score = compute_spatial_chaos(pred.values, spatial_coords)    
    asw_score = compute_asw_spatial(adata, "predict_cluster", coords=spatial_coords)
    scs_score = compute_spatial_coherence_score(labels=pred.values, coords=spatial_coords, n_neighbors=10,n_permutations=100, random_state=0)
    
    print(pas_score, chaos_score, asw_score, scs_score)
    print(f"Done. Continue..")

    # == Save summary metrics ==
    print(f"Saving results for ...")
    results_summary.append({
        'file_name': f"{file_predicted_clustering}_{method}_perf_and_spatial_summary.csv", 
        'config': in_file,
        'ARI': round(ari, 4),
        'NMI': round(nmi, 4),
        'FMI': round(fmi, 4),
        'Purity': round(purity, 4),
        'Homogeneity': round(homog, 4),
        'Completeness': round(comp, 4),
        'Vmeasure': round(vmeas, 4),
        'PAS': round(pas_score, 4),
        'CHAOS': round(chaos_score, 4),
        'ASW_fix': round(asw_score, 4),
        'SCS':round(scs_score, 4),
        'n_clusters_pred': n_clusters_pred,
        'n_clusters_gt': n_clusters_gt,
        'n_cells': len(adata)
    })

    summary_df = pd.DataFrame(results_summary)
    summary_df.to_csv(csv_path, index=False)

    return summary_df

import argparse
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_file', required=True)
    parser.add_argument('--out_path', type=str, required=True)
    parser.add_argument('--meta_dir', type=str, required=True)
    parser.add_argument('--file_predicted_clustering', required=True)
    parser.add_argument('--method', type=str, default='Banksy')
    parser.add_argument('--annotation_key', required=True)
    parser.add_argument('--spatial_x', required=True)
    parser.add_argument('--spatial_y', required=True)

    args = parser.parse_args()

    main(
        in_file=args.in_file,
        out_path=Path(args.out_path),
        meta_dir=Path(args.meta_dir),
        file_predicted_clustering=args.file_predicted_clustering,
        method=args.method,
        annotation_key=args.annotation_key,
        spatial_x=args.spatial_x,
        spatial_y=args.spatial_y
    )
