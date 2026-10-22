library(argparse)
library(BASS)
library(SingleCellExperiment)
library(Matrix)
library(scater)
library(dplyr)
library(anndata)
library(pryr)
library(fs)

# == Define helpers ==
params_to_string <- function(k_geom, pca_dims, n_hvg, init_method) {
  pca_str <- paste(pca_dims, collapse = "-")
  return(paste0("hvg", n_hvg,"_method_", init_method, "_k", k_geom, "_pca", pca_str))
}

# == Define main ==
main <- function(in_file, annotation_key, x_spatial, y_spatial, init_method = "kmeans", k_geom, n_hvg=3000, pca_dims=20, iterations = 10){
    # Track time and memory
    gc(reset = TRUE)
    start_time <- Sys.time()
    mem_start <- sum(gc()[, 2])
    message(start_time)

    # == Paths ==
    in_path <- file.path( "/faststorage/project/limitedAccessBladder/NOBACKUP/BIRC_thesis_students_spring26/marinarb/outputs/normalized_log1p_hvgs/", paste0(in_file, ".h5ad"))
   
    out_path <- file.path( "/faststorage/project/limitedAccessBladder/NOBACKUP/BIRC_thesis_students_spring26/marinarb/outputs/BASS/")

    # == Load ==
    # the data imported is already normalized by depth, log1p transformed and feature selected

    message("\n Loading anndata object...")
    adata <- read_h5ad(in_path)
    cat("Shape:", nrow(adata$X), "x", ncol(adata$X), "\n")
    
    message("\n Loaded:", in_path)

    # Conversion to SCE object
    message("\n Conversion to SCE...")
    
    X_matrix <- as.matrix(adata$X)
    class(X_matrix)
    message("Dimension: ", paste(dim(X_matrix), collapse = " x "))
    
    sce <- SingleCellExperiment(
        assays = list(counts = t(X_matrix)),
        colData = DataFrame(adata$obs),
        rowData = DataFrame(adata$var))
    
    # Sanity: matching dimensions
    stopifnot(nrow(assay(sce, "counts")) == nrow(adata$var))  # genes
    stopifnot(ncol(assay(sce, "counts")) == nrow(adata$obs))  # cells
    
    message("\n Overview of SCE:")
    print(sce)
    
    # Define spatial coordinates
    if (x_spatial %in% colnames(adata$obs) && y_spatial %in% colnames(adata$obs)) {
        spatial_coords <- as.matrix(adata$obs[, c(x_spatial, y_spatial)])
        
        reducedDim(sce, "spatial_bass") <- spatial_coords
        
        colData(sce)$x_bass <- spatial_coords[, 1]
        colData(sce)$y_bass <- spatial_coords[, 2]
        
        message("Loaded spatial coordinates: x_bass, y_bass, spatial_bass")
    } else {
        stop("Spatial coordinates not found: ", x_spatial, "/", y_spatial)  # más estricto
    }
    
    coord_keys <- c("x_bass", "y_bass", "spatial_bass")
    message("Coordinate keys: ", paste(coord_keys, collapse = ", "))
    
    # Ground truth
    ground_truth <- colData(sce)[[annotation_key]]  
    n_gt_clusters <- length(unique(ground_truth[!is.na(ground_truth)]))
    message("Ground truth clusters (", annotation_key, "): ", n_gt_clusters)
    
    message("\nSolved. Continue...")
    
    # -- Extract HVGs already selected with Scanpy --
    message("Extract HVGs...")
    
    if ("highly_variable" %in% colnames(rowData(sce))) {
        hvgs <- rowData(sce)$highly_variable
    } else if ("highly_variable_genes" %in% colnames(rowData(sce))) {
        hvgs <- rowData(sce)$highly_variable_genes
    } else {
        stop("No HVG column found in rowData(sce)")
    }
    
    hvg_indices <- which(hvgs)
    if(length(hvg_indices) == 0) {
        stop("No HVGs found in rowData(sce)")
    }
    if(length(hvg_indices) > n_hvg) {
        var_hvgs <- apply(assay(sce, "counts")[hvg_indices, ], 1, var)
        hvg_indices <- hvg_indices[order(var_hvgs, decreasing = TRUE)[1:n_hvg]]
    }
        
    message("\nSolved. Continue...")

    # == RUNNING BASS ==
    # == Create BASS object ==
    counts <- assays(sce)$counts
    counts <- counts[hvg_indices, ] # force Scanpy HVGs
    xym <- as.matrix(colData(sce)[, c("x_bass", "y_bass")]) 
    cntm <- list()
    cntm[[1]] <- counts
    xym_list <- list()
    xym_list[[1]] <- xym 
    
    C <- n_gt_clusters  # celltype clusters
    R <- 6   # spatial domains

    # Store results in meta
    cell_candidates <- c("cell_key", "cell_id", "cell")
    cell_col <- NULL
    
    for(candidate in cell_candidates) {
        if(candidate %in% colnames(colData(sce))) {
            cell_col <- candidate
            break
        }
    }
    
    if(!is.null(cell_col)) {
        cell_id <- colData(sce)[[cell_col]]
        message("Using cell ID column: ", cell_col)
    } else {
        cell_id <- colnames(sce)
        message("No cell ID column found, using colnames(sce)")
    }
    
    meta <- data.frame(cell_key = as.character(cell_id))
    
    set.seed(0)
    
    BASS <- BASS:::createBASSObject(
        cntm,
        xym_list,
        C,
        R,
        beta_method = "SW",
        k = k_geom, # explicitly define k-NN neighbours to use
        init_method = init_method,
        nsample = iterations # MCMC iterations for cluster refinement - reduce for big datasets
    )

    # == Run BASS pre-processing ==
    BASS <- BASS:::BASS.preprocess(
      BASS,
      doLogNormalize = FALSE,
      doPCA = TRUE,
      scaleFeature = FALSE,
      doBatchCorrect = FALSE,
      nHVG = n_hvg,
      nPC = pca_dims
    )
    
    message("BASS hyperparameters:")
    listAllHyper(BASS)

    # == Run BASS algorithm ==
    message("Starting BASS...")
    
    BASS <- BASS.run(BASS)
    table(BASS@results$c[[1]])
    colData(sce)$BASS <- BASS@results$c[[1]]
    
    message("BASS clustering completed")
    
    # == Run BASS post-processing ==
    message("Starting BASS postprocessing...")
    
    BASS <- BASS.postprocess(BASS)
    
    message("Postprocessing completed")
    
    # == Extract results ==
    zlabels <- BASS@results$z # predicted spatial domain labels
    clabels <- BASS@results$c # predicted cell type labels
    
    meta$"predict_cluster" <- clabels[[1]]   
     
    # Save
    param_str <- params_to_string(k_geom, pca_dims, n_hvg, init_method)
    out_file <- paste0(in_file, "_", param_str,"_BASS_predclusters" ,".csv")
    
    csv_path <- file.path(out_path, out_file)
    write.csv(meta, csv_path, row.names = FALSE)
    message(paste0("Results for ", in_file, " saved in ", csv_path))

    # Tack
    end_time <- Sys.time()
    message(end_time)
    runtime <- end_time - start_time
    gc()
    mem_end <- sum(gc()[, 2])
    used_memory <- mem_end - mem_start
    message(paste0("Execution time:", runtime))
    message("Peak memory used: ", round(used_memory, 1), " MB") 
}

parser <- ArgumentParser(description = "Run BASS delltype spatial clustering")
parser$add_argument('--in_file', required = TRUE)
parser$add_argument('--annotation_key', type = "character", required = TRUE)
parser$add_argument('--x_spatial', type = "character", required = TRUE)
parser$add_argument('--y_spatial', type = "character", required = TRUE)
parser$add_argument('--init_method', type = "character", default = "kmeans")
parser$add_argument('--k_geom', type = "integer", default = 15)
parser$add_argument('--n_hvg', type = "integer", default = 3000)
parser$add_argument('--pca_dims', type = "integer", default = 20)
parser$add_argument('--iterations', type = "integer", default = 10)

args <- parser$parse_args()
main(args$in_file, args$annotation_key, args$x_spatial, args$y_spatial, args$init_method, args$k_geom, args$n_hvg, args$pca_dims, args$iterations)
