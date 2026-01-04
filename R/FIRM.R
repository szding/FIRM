#' @useDynLib FIRM, .registration = TRUE
NULL

#' Flexible Integration of Single-Cell RNA-Seq Data
#'
#' Performs unsupervised integration of two single-cell RNA-seq datasets
#' by searching for the optimal clustering resolution
#' pair that maximises mutual-nearest-neighbour (MNN) mixing in the
#' combined PCA space.  The final integrated expression matrix is returned
#' after batch-effect correction.
#'
#' @param SS2 \code{Seurat} object for reference dataset (e.g. Smart-seq2).
#' @param tenx \code{Seurat} object for query dataset (e.g. 10x Genomics).
#' @param hvg1,hvg2 Character vectors giving the high-variable gene names
#'   selected in SS2 and tenx, respectively.
#' @param dims Integer scalar, number of principal components to use
#'   during integration.
#' @param all_genes Logical scalar.  If \code{TRUE} the integration is
#'   carried out on the union of all genes; otherwise only on the
#'   intersected HVGs (default).
#' @param res_seq_SS2,res_seq_tenx Numeric vectors of clustering
#'   resolutions to be screened for the reference and query dataset,
#'   respectively.  Defaults to \code{seq(0.1, 2, 0.1)}.
#' @param coreNum Integer scalar, number of CPU cores used for parallel
#'   screening.
#' @param verbose Logical scalar.  If \code{TRUE} a list containing the
#'   integrated matrix and quality metrics is returned; otherwise only
#'   the integrated expression matrix is returned.
#'
#' @return
#' By default (\code{verbose = FALSE}) a single \code{matrix} of
#' batch-corrected, scaled expression values with genes as rows and
#' combined cells as columns.
#'
#' If \code{verbose = TRUE} a named \code{list} is returned:
#' \describe{
#'   \item{integrated}{As above, the corrected expression matrix.}
#'   \item{Metric_PCA}{Mean MNN mixing score of the naive PCA
#'     (no correction).}
#'   \item{Metric_FIRM}{Mean MNN mixing score of the FIRM-corrected
#'     embedding (matrix when multiple resolution pairs were screened).}
#' }
#'
#' @details
#' The algorithm performs the following steps:
#' \enumerate{
#'   \item PCA on each dataset using the intersected HVGs.
#'   \item SNN graph construction (via \code{Seurat::FindNeighbors}).
#'   \item Screening of clustering resolution pairs
#'     (\code{res_seq_SS2} × \code{res_seq_tenx}) to maximise
#'     mutual-nearest-neighbour mixing in the joint PCA space.
#'   \item Batch-effect correction with \code{FIRM_res*} functions.
#'   \item Final integrated expression matrix is scaled and returned.
#' }
#'
#' Quality control: the integrated embedding is compared with the naive
#' PCA; if correction does not improve mixing the latter is returned.
#'
#' @seealso
#' \code{\link{prep_data}} for Data preprocessing.
#'
#' @references
#' Ming, J., Lin, Z., Zhao, J., Wan, X., Ezran, C., Liu, S., ... & TTM Consortium. (2022). FIRM: Flexible integration of single-cell RNA-sequencing data for large-scale multi-tissue cell atlas datasets. \emph{Briefings in bioinformatics}, 23(5).
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' data("ExampleData")
#' SS2  <- ExampleData$SS2
#' tenx <- ExampleData$tenx
#' hvg1 <- VariableFeatures(SS2)
#' hvg2 <- VariableFeatures(tenx)
#'
#' integrated <- FIRM(SS2, tenx, hvg1, hvg2, dims = 30)
#' dim(integrated)
#'
#' # with full diagnostics
#' res <- FIRM(SS2, tenx, hvg1, hvg2, dims = 30, verbose = TRUE)
#' str(res)
#' }
#'
#' @importFrom RANN nn2
#' @importFrom Seurat CreateSeuratObject NormalizeData FindVariableFeatures
#'   ScaleData VariableFeatures RunPCA FindNeighbors FindClusters GetAssayData
#'   SetAssayData Embeddings
#' @export
FIRM <- function(SS2, tenx, hvg1, hvg2, dims, all_genes = FALSE, res_seq_SS2 = seq(0.1, 2, 0.1),
                 res_seq_tenx = seq(0.1, 2, 0.1),
                 coreNum = 1, verbose = FALSE){

  set.seed(0)
  hvg <- intersect(hvg1, hvg2)
  gene_all <- union(rownames(SS2), rownames(tenx))

  SS2 <- RunPCA(SS2, features = hvg, npcs = dims, verbose = FALSE)
  tenx <- RunPCA(tenx, features = hvg, npcs = dims, verbose = FALSE)
  gc()

  SS2 <- FindNeighbors(SS2, reduction = "pca", dims = 1:dims)
  SS2_snn <- SS2[["RNA_snn"]]
  tenx <- FindNeighbors(tenx, reduction = "pca", dims = 1:dims)
  tenx_snn <- tenx[["RNA_snn"]]

  SS2 <- GetAssayData(SS2, assay = "RNA", layer = "scale.data")
  tenx <- GetAssayData(tenx, assay = "RNA", layer = "scale.data")


  gc()

  res_SS2 <- NULL
  res_tenx <- NULL
  for (i in 1:length(res_seq_SS2)){
    if (is.null(tryCatch(SS2_FindClusters <- FindClusters(SS2_snn, resolution = res_seq_SS2[i], verbose = FALSE)[, 1], error = function(e){}))){
      break
    }

    SS2_FindClusters <- FindClusters(SS2_snn, resolution = res_seq_SS2[i], verbose = FALSE)[, 1]
    num_cluster_SS2 <- length(unique(SS2_FindClusters))

    if (num_cluster_SS2 > 1){
      if (is.null(res_SS2)){
        res_SS2 <- c(res_SS2, res_seq_SS2[i])
        num_cluster_SS2_old <- length(unique(SS2_FindClusters))
      }
      else if (num_cluster_SS2 > num_cluster_SS2_old){
        res_SS2 <- c(res_SS2, res_seq_SS2[i])
      }
    }
    num_cluster_SS2_old <- num_cluster_SS2
  }

  if (all_genes == FALSE){
    integrated_PCA <- cbind(SS2[hvg, ], tenx[hvg, ])
  } else {
    integrated_PCA <- matrix(0, length(gene_all), ncol(SS2) + ncol(tenx))
    rownames(integrated_PCA) <- gene_all
    colnames(integrated_PCA) <- c(colnames(SS2), colnames(tenx))
    integrated_PCA[rownames(SS2), 1:ncol(SS2)] <- SS2
    integrated_PCA[rownames(tenx), (ncol(SS2)+1):(ncol(SS2)+ncol(tenx))] <- tenx
  }

  integrated_PCA_obj <- CreateSeuratObject(counts = integrated_PCA)
  integrated_PCA_obj <- NormalizeData(integrated_PCA_obj)
  integrated_PCA_obj <- SetAssayData(object = integrated_PCA_obj, assay = "RNA", layer = "data", new.data = integrated_PCA)
  integrated_PCA_obj <- ScaleData(integrated_PCA_obj, do.center = FALSE, verbose = FALSE)
  integrated_PCA <- GetAssayData(integrated_PCA_obj, assay = "RNA", layer = "scale.data")


  dataset_list_PCA <- c(rep(1, ncol(SS2)), rep(2, ncol(tenx)))

  if (length(res_SS2) == 0){
    if (verbose == TRUE){
      return(list(integrated = integrated_PCA))
    }
    else{
      return(integrated_PCA)
    }
  }

  for (i in 1:length(res_seq_tenx)){
    if (is.null(tryCatch(tenx_FindClusters <- FindClusters(tenx_snn, resolution = res_seq_tenx[i], verbose = FALSE)[, 1], error = function(e){}))){
      break
    }

    tenx_FindClusters <- FindClusters(tenx_snn, resolution = res_seq_tenx[i], verbose = FALSE)[, 1]
    num_cluster_tenx <- length(unique(tenx_FindClusters))

    if (num_cluster_tenx > 1){
      if (is.null(res_tenx)){
        res_tenx <- c(res_tenx, res_seq_tenx[i])
        num_cluster_tenx_old <- length(unique(tenx_FindClusters))
      }
      else if (num_cluster_tenx > num_cluster_tenx_old){
        res_tenx <- c(res_tenx, res_seq_tenx[i])
      }
    }
    num_cluster_tenx_old <- num_cluster_tenx
  }

  res_SS2 <- res_SS2[length(res_SS2):1]
  res_tenx <- res_tenx[length(res_tenx):1]

  SS2_FindClusters <- FindClusters(SS2_snn, resolution = res_SS2, verbose = FALSE)
  tenx_FindClusters <- FindClusters(tenx_snn, resolution = res_tenx, verbose = FALSE)

  rm(SS2_snn)
  rm(tenx_snn)
  gc()

  SS2_FindClusters <- matrix(as.numeric(as.matrix(SS2_FindClusters)), nrow(SS2_FindClusters), ncol(SS2_FindClusters))
  tenx_FindClusters <- matrix(as.numeric(as.matrix(tenx_FindClusters)), nrow(tenx_FindClusters), ncol(tenx_FindClusters))

  if (all_genes == FALSE){
    if (ncol(SS2) < ncol(tenx)){
      Dataset1 <- SS2[hvg, ]
      Dataset2 <- tenx[hvg, ]
      FindClusters1 <- SS2_FindClusters
      FindClusters2 <- tenx_FindClusters
      res1 <- res_SS2
      res2 <- res_tenx
    } else {
      Dataset1 <- tenx[hvg, ]
      Dataset2 <- SS2[hvg, ]
      FindClusters1 <- tenx_FindClusters
      FindClusters2 <- SS2_FindClusters
      res1 <- res_tenx
      res2 <- res_SS2
    }
  } else {
    gene_all_num <- length(gene_all)

    tmp <- seq(1, gene_all_num, 1)
    names(tmp) <- gene_all
    gene_all_ind_SS2 <- as.numeric(tmp[rownames(SS2)[which(rownames(SS2) %in% gene_all)]])
    tmp <- seq(1, nrow(SS2), 1)
    names(tmp) <- rownames(SS2)
    hvg_ind_SS2 <- as.numeric(tmp[hvg])

    tmp <- seq(1, gene_all_num, 1)
    names(tmp) <- gene_all
    gene_all_ind_tenx <- as.numeric(tmp[rownames(tenx)[which(rownames(tenx) %in% gene_all)]])
    tmp <- seq(1, nrow(tenx), 1)
    names(tmp) <- rownames(tenx)
    hvg_ind_tenx <- as.numeric(tmp[hvg])

    gene_all_hvg <- which(gene_all %in% hvg)

    if (ncol(SS2) < ncol(tenx)){
      Dataset1 <- SS2
      Dataset2 <- tenx
      hvg_ind1 <- hvg_ind_SS2
      FindClusters1 <- SS2_FindClusters
      hvg_ind2 <- hvg_ind_tenx
      FindClusters2 <- tenx_FindClusters
      gene_all_ind1 <- gene_all_ind_SS2
      gene_all_ind2 <- gene_all_ind_tenx
      res1 <- res_SS2
      res2 <- res_tenx
    } else {
      Dataset1 <- tenx
      Dataset2 <- SS2
      hvg_ind1 <- hvg_ind_tenx
      FindClusters1 <- tenx_FindClusters
      hvg_ind2 <- hvg_ind_SS2
      FindClusters2 <- SS2_FindClusters
      gene_all_ind1 <- gene_all_ind_tenx
      gene_all_ind2 <- gene_all_ind_SS2
      res1 <- res_tenx
      res2 <- res_SS2
    }
  }

  dataset_list <- c(rep(1, ncol(Dataset1)), rep(2, ncol(Dataset2)))
  quantile_default <- 0.75
  max.k <- 300
  rept <- 50

  if ((length(res1)*length(res2)) == 1){
    if (all_genes == FALSE){
      result <- FIRM_res_hvg(Dataset1, FindClusters1, Dataset2, FindClusters2,
                         dims, quantile_default = quantile_default,rept_ds = rept)
    } else {
      result <- FIRM_res(Dataset1, hvg_ind1, FindClusters1,
                         Dataset2, hvg_ind2, FindClusters2,
                         dims, gene_all_num, gene_all_hvg,
                         gene_all_ind1, gene_all_ind2,
                         quantile_default = quantile_default,rept_ds = rept)
    }

    integrated_PCA_obj <- RunPCA(integrated_PCA_obj, features = hvg, npcs = dims, verbose = FALSE)
    Metric_PCA <- mean(Mixing_Metric(Embeddings(integrated_PCA_obj, "pca"), dataset_list_PCA, max.k = max.k))
    rm(integrated_PCA_obj)

    if (all(result$integrated == 0)){
      if (verbose == TRUE){
        return(list(integrated = integrated_PCA, Metric_PCA = Metric_PCA))
      }
      else{
        return(integrated_PCA)
      }
    }

    integrated_FIRM <- result$integrated
    if (all_genes == FALSE){
      rownames(integrated_FIRM) <- hvg
    } else {
      rownames(integrated_FIRM) <- gene_all
    }
    colnames(integrated_FIRM) <- c(colnames(Dataset1), colnames(Dataset2))

    integrated_FIRM_obj <- CreateSeuratObject(counts = integrated_FIRM)
    integrated_FIRM_obj <- NormalizeData(integrated_FIRM_obj)
    integrated_FIRM_obj <- SetAssayData(object = integrated_FIRM_obj, assay = "RNA", layer = "data", new.data = integrated_FIRM)
    integrated_FIRM_obj <- ScaleData(integrated_FIRM_obj, do.center = FALSE, verbose = FALSE)
    integrated_FIRM <- GetAssayData(integrated_FIRM_obj, assay = "RNA", layer = "scale.data")
    integrated_FIRM_obj <- RunPCA(integrated_FIRM[hvg, ], features = hvg, npcs = dims, verbose = FALSE)

    pca_embeddings <- Embeddings(integrated_FIRM_obj, reduction = "pca")
    rm(integrated_FIRM_obj)
    Metric_FIRM <- mean(Mixing_Metric(pca_embeddings, dataset_list, max.k = max.k))

    gc()

    if(min(Metric_FIRM) >= Metric_PCA){
      if (verbose == TRUE){
        return(list(integrated = integrated_PCA, Metric_PCA = Metric_PCA, Metric_FIRM = Metric_FIRM))
      }
      else{
        return(integrated_PCA)
      }
    } else {
      if (verbose == TRUE){
        return(list(integrated = integrated_FIRM, Metric_FIRM = Metric_FIRM, Metric_PCA = Metric_PCA))
      }
      else{
        return(integrated_FIRM)
      }
    }
  }

  result <- FIRM_res_all(Dataset1[hvg, ], FindClusters1, Dataset2[hvg, ], FindClusters2,
                         dims, quantile_default = quantile_default,rept_ds = rept, coreNum = coreNum)

  Metric_PCA <- mean(Mixing_Metric(result$Embedding_PCA, dataset_list, max.k = max.k))

  Metric_FIRM <- matrix(0, length(res1), length(res2))
  rownames(Metric_FIRM) <- res1
  colnames(Metric_FIRM) <- res2

  idx <- 0;
  for (i in 1:length(res1)){
    for (j in 1:length(res2)){
      idx <- idx + 1
      if (all(result$Embedding_FIRM[, , idx] == 0)){
        Metric_FIRM[i, j] <- Metric_PCA
      } else {
        Metric_FIRM[i, j] <- mean(Mixing_Metric(result$Embedding_FIRM[, , idx], dataset_list, max.k = max.k))
      }
    }
  }

  if(min(Metric_FIRM) >= Metric_PCA){
    if (verbose == TRUE){
      return(list(integrated = integrated_PCA, Metric_PCA = Metric_PCA, Metric_FIRM = Metric_FIRM))
    }
    else{
      return(integrated_PCA)
    }
  } else{
    i <- which.min(Metric_FIRM) %% length(res1)
    if (i == 0){
      i <- length(res1)
    }
    j <- ceiling(which.min(Metric_FIRM)/length(res1))

    if (all_genes == FALSE){
      result <- FIRM_res_hvg(Dataset1, FindClusters1[, i], Dataset2, FindClusters2[, j],
                         dims, quantile_default = quantile_default,rept_ds = rept)
    } else {
      result <- FIRM_res(Dataset1, hvg_ind1, FindClusters1[, i],
                         Dataset2, hvg_ind2, FindClusters2[, j],
                         dims, gene_all_num, gene_all_hvg,
                         gene_all_ind1, gene_all_ind2,
                         quantile_default = quantile_default,rept_ds = rept)
    }

    if (all(result$integrated == 0)){
      if (verbose == TRUE){
        return(list(integrated = integrated_PCA, Metric_PCA = Metric_PCA, Metric_FIRM = Metric_FIRM))
      }
      else{
        return(integrated_PCA)
      }
    }

    integrated_FIRM <- result$integrated
    if (all_genes == FALSE) {
      rownames(integrated_FIRM) <- hvg
    } else {
      rownames(integrated_FIRM) <- gene_all
    }
    colnames(integrated_FIRM) <- c(colnames(Dataset1), colnames(Dataset2))

    integrated_FIRM_obj <- CreateSeuratObject(counts = integrated_FIRM)
    integrated_FIRM_obj <- NormalizeData(integrated_FIRM_obj)
    integrated_FIRM_obj <- SetAssayData(object = integrated_FIRM_obj, assay = "RNA", layer = "data", new.data = integrated_FIRM)
    integrated_FIRM_obj <- ScaleData(integrated_FIRM_obj, do.center = FALSE, verbose = FALSE)
    integrated_FIRM <- GetAssayData(integrated_FIRM_obj, assay = "RNA", layer = "scale.data")
    rm(integrated_FIRM_obj)

    gc()

    if (verbose == TRUE){
      return(list(integrated = integrated_FIRM, Metric_FIRM = Metric_FIRM, Metric_PCA = Metric_PCA))
    }
    else{
      return(integrated_FIRM)
    }
  }
}

