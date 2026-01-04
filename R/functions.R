#' Data preprocessing for FIRM integration
#'
#' Performs a standard Seurat workflow: normalization, scaling and
#' selection of the top 4 000 highly-variable genes (HVGs).
#'
#' @param counts          Raw count matrix (genes × cells) or \code{dgCMatrix}.
#' @param file_path       Optional directory where results should be **saved**.
#'                        If \code{NULL} (default) nothing is written to disk.
#' @param file_name       Optional basename for saved \code{.RData} files.
#'                        When given, two files are created:
#'                        \code{<file_name>.RData} (Seurat object) and
#'                        \code{<file_name>_hvg.RData} (character vector of HVGs).
#' @param hvg_genes       Target number of genes to return (default 4 000).
#'
#' @return  A named list with elements
#'          \describe{
#'            \item{Dataset}{Seurat object after normalization, scaling and
#'                           feature selection.}
#'            \item{hvg}{Character vector of the 4 000 most variable genes.}
#'          }
#'
#' @author  Jingsi Ming
#'
#' @seealso  \code{\link{FIRM}} for the integration step.
#'
#' @examples
#' \dontrun{
#' data("ExampleData")
#' prep_SS2  <- prep_data(ExampleData$SS2, hvg_genes = 1000)
#' Dataset1  <- prep_SS2$Dataset
#' hvg1      <- prep_SS2$hvg
#'
#' prep_tenx <- prep_data(ExampleData$tenx, hvg_genes = 1000)
#' Dataset2  <- prep_tenx$Dataset
#' hvg2      <- prep_tenx$hvg
#'
#' # save results to disk
#' prep_data(ExampleData$tenx,
#'           file_path = tempdir(),
#'           file_name = "tenx_processed")
#' }
#'
#' @importFrom Seurat CreateSeuratObject NormalizeData FindVariableFeatures
#'   ScaleData VariableFeatures
#' @export
prep_data <- function(counts, file_path = NULL, file_name = NULL, hvg_genes = 4000) {
  Dataset <- CreateSeuratObject(counts = counts)
  Dataset <- NormalizeData(Dataset, layer = "counts", verbose = FALSE)
  Dataset <- FindVariableFeatures(Dataset, nfeatures = hvg_genes,
                                  layer = "data", verbose = FALSE)
  Dataset <- ScaleData(Dataset, features = rownames(Dataset),
                       layer = "data", do.center = FALSE, verbose = FALSE)

  hvg <- VariableFeatures(Dataset)

  return(list(Dataset = Dataset, hvg = hvg))
}


#' Select consensus genes across multiple HVG lists
#'
#' Rank-based selection of genes that appear in multiple HVG lists,
#' prioritising those present in all datasets.
#'
#' @param hvg_list  List of character vectors, each containing HVGs from one dataset.
#' @param gene_all  Optional character vector of all candidate genes (intersection filter).
#' @param num       Target number of genes to return (default 4 000).
#'
#' @return Character vector of selected genes.
#'
#' @importFrom stats median
#' @export
SelectGene <- function(hvg_list, gene_all = NULL, num = 4000){
  K <- length(hvg_list)

  hvg_union <- NULL
  for (k in 1:K){
    hvg_union <- union(hvg_union, hvg_list[[k]])
    if (!is.null(gene_all)){
      hvg_union <- intersect(hvg_union, gene_all)
    }
  }

  if (length(hvg_union) < num){
    return(hvg_union)
  }

  hvg_union_rank <- matrix(0, length(hvg_union), K)
  rownames(hvg_union_rank) <- hvg_union
  for (i in 1:length(hvg_union)){
    for (k in 1:K){
      if (hvg_union[i] %in% hvg_list[[k]]) {
        hvg_union_rank[i, k] <- which(hvg_list[[k]] == hvg_union[i])
      }
    }
  }

  hvg <- hvg_union[which(rowSums(hvg_union_rank != 0) == K)]

  if (length(hvg) >= num){
    return(hvg)
  } else{
    for (k in (K-1):1){
      hvg_add <- hvg_union[which(rowSums(hvg_union_rank != 0) == k)]
      if (length(hvg_add) + length(hvg) >= num){
        hvg_add <- names(sort(apply(hvg_union_rank[hvg_add, ], 1, function(x) stats::median(x[x != 0]))))[1:(num-length(hvg))]
        hvg <- c(hvg, hvg_add)
        return(hvg)
      } else {
        hvg <- c(hvg, hvg_add)
      }
    }
  }
}

#' Combine or select top HVGs from saved prep_data outputs
#'
#' Load HVG lists saved by \code{prep_data} and either take their union
#' or rank-based top consensus.
#'
#' @param file_names Character vector of basenames (without "_hvg.RData").
#' @param file_path  Directory containing the \code{*_hvg.RData} files.
#' @param method     "all" for union, "top" for ranked consensus (see \code{SelectGene}).
#' @param hvg_genes  Target number of genes to return (default 4 000).
#'
#' @return Character vector of selected genes.
#'
#' @seealso \code{\link{prep_data}}
#' @export
Select_hvg <- function(file_names, file_path, method, hvg_genes = 4000){
  hvg <- NULL
  if (method == "all"){
    hvg_all <- NULL
    for (i in 1:length(file_names)){
      load(paste(file_path, "/", file_names[i], "_hvg.RData", sep = ""))
      hvg_all <- union(hvg_all, hvg)
    }
    return(hvg_all)
  }

  if (method == "top"){
    hvg_all <- NULL
    for (i in 1:length(file_names)){
      load(paste(file_path, "/", file_names[i], "_hvg.RData", sep = ""))
      hvg_all <- c(hvg_all, list(hvg))
    }
    hvg_top <- SelectGene(hvg_all, num = hvg_genes)
    return(hvg_top)
  }
}
