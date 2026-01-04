#' Example single-cell datasets
#'
#' A list of Seurat objects and corresponding metadata from two 10x Genomics
#' peripheral blood samples (SS2 and tenx).
#'
#' @format A list of 4 elements:
#' \describe{
#'   \item{SS2}{Seurat object from Smart-seq2 platform}
#'   \item{tenx}{Seurat object from 10x Chromium platform}
#'   \item{meta_SS2}{data.frame, cell-level metadata for SS2}
#'   \item{meta_tenx}{data.frame, cell-level metadata for tenx}
#' }
#' @examples
#' \dontrun{
#' data("ExampleData")
#' names(ExampleData)
#' head(ExampleData$meta_SS2)
#' }
#' @keywords datasets
"ExampleData"


#' k-NN Mixing Metric
#'
#' Evaluate the degree of mixing between two or more datasets in a low-dimensional embedding.
#'
#' @param embedding       Numeric matrix (n × d) where each row is the coordinates of one
#'                        sample in the low-dimensional space (PCA, UMAP, t-SNE, …).
#' @param dataset_list    Factor or character vector of length n indicating which dataset
#'                        each sample comes from.
#' @param k               Positive integer.  The rank of the within-dataset neighbour to look for
#'                        (default 5).
#' @param max.k           Positive integer.  Total number of nearest neighbours to compute
#'                        (default 300).
#' @return                Named numeric vector with one entry per unique level of
#'                        `dataset_list`.  The entry is the median position (among
#'                        1 … max.k) of the k-th within-dataset neighbour across all
#'                        samples that belong to that dataset.  Values are clamped to
#'                        max.k when fewer than k within-dataset neighbours exist.
#'
#' @importFrom stats median
#' @rdname Mixing_Metric
#' @export
Mixing_Metric <- function(embedding, dataset_list, k = 5, max.k = 300) {
  set.seed(0)
  nn <- RANN::nn2(embedding, k = max.k)$nn.idx[, -1]

  mixing <- sapply(
    X = 1:nrow(x = nn),
    FUN = function(x) {
      sapply(X = unique(dataset_list), FUN = function(y) {
        which(x = dataset_list[nn[x, ]] == y)[k]
      })
    }
  )

  mixing[is.na(x = mixing)] <- max.k
  mixing <- apply(
    X = mixing,
    MARGIN = 2,
    FUN = stats::median
  )

  return(mixing)
}


#' Local structure preservation metric
#'
#' Compute the neighbourhood overlap between the original PCA space of each
#' dataset and the integrated PCA space, averaged over the top \code{neighbors}
#' nearest neighbours.  High overlap indicates better local structure
#' preservation after integration.
#'
#' @param SS2,tenx  \link[Seurat]{Seurat} objects containing the **unintegrated**
#'   PCA reduction (slot \code{reductions$pca}).
#' @param integrated  \link[Seurat]{Seurat} object that contains the **joint**
#'   PCA reduction and a meta-data column \code{dataset}.
#' @param dims  Number of principal components to use (default 30).
#' @param neighbors  Size of the neighbourhood to evaluate (default 20).
#' @param emb_key  Name of the reduction to use (default \code{"pca"}).
#' @param batch_key  Name of the meta-data column that indicates dataset/batch
#'   origin (default \code{"dataset"}).
#'
#' @return
#' A named list with two numeric vectors:
#' \item{SS2}{Per-cell overlap scores for the SS2 dataset.}
#' \item{tenx}{Per-cell overlap scores for the 10X dataset.}
#'
#' @rdname Local_Struct
#' @export
Local_Struct <- function(SS2, tenx, integrated, dims = 30, emb_key = "pca",
                         batch_key = "dataset", neighbors = 20) {
  nn.SS2 <- RANN::nn2(Embeddings(SS2, reduction = emb_key)[, 1:dims], k = neighbors)$nn.idx
  nn.tenx <- RANN::nn2(Embeddings(tenx, reduction = emb_key)[, 1:dims], k = neighbors)$nn.idx

  ss2_cells <- which(integrated[[batch_key]] == unique(SS2[[batch_key]]))
  tenx_cells <- which(integrated[[batch_key]] == unique(tenx[[batch_key]]))
  nn.corrected.SS2 <- RANN::nn2(Embeddings(integrated, reduction = emb_key)[ss2_cells, 1:dims], k = neighbors)$nn.idx
  nn.corrected.tenx <- RANN::nn2(Embeddings(integrated, reduction = emb_key)[tenx_cells, 1:dims], k = neighbors)$nn.idx

  local.struct.SS2 <- sapply(X = 1:nrow(x = nn.SS2), FUN = function(x) {
    length(x = intersect(x = nn.SS2[x, ], y = nn.corrected.SS2[x, ])) / neighbors})
  local.struct.tenx <- sapply(X = 1:nrow(x = nn.tenx), FUN = function(x) {
    length(x = intersect(x = nn.tenx[x, ], y = nn.corrected.tenx[x, ])) / neighbors})

  local.struct <- list(SS2 = local.struct.SS2, tenx = local.struct.tenx)
  return(local.struct)
}

#' Label transfer via k-NN (hard assignment)
#'
#' For every query cell, find its \code{k} nearest reference cells in the
#' integrated PCA space and assign the most frequent label among them.
#'
#' @param integrated  \link[Seurat]{Seurat} object with integrated PCA and
#'   meta-data columns \code{dataset} and \code{annotation}.
#' @param ref  Character string identifying the reference dataset (default
#'   \code{"10X"}).
#' @param query  Character string identifying the query dataset (default
#'   \code{"SS2"}).
#' @param label_key  Name of the meta-data column that indicates celltype
#' @param emb_key  Name of the reduction to use (default \code{"pca"}).
#' @param batch_key  Name of the meta-data column that indicates dataset/batch
#' @param k  Number of nearest neighbours to vote (default 10).
#'
#' @return
#' Named character vector: predicted label for each query cell (names are
#' cell barcodes).
#'
#' @rdname label_trans
#' @export
label_trans <- function(integrated, ref = "10X", query = "SS2", label_key = "annotation",
                        emb_key = "pca", batch_key = "dataset", k = 10){
  set.seed(0)

  anno_tenx <- integrated[[label_key]][which(integrated[[batch_key]] == ref)]

  ref_embeddings <- Embeddings(integrated, reduction = emb_key)[which(integrated[[batch_key]] == ref), ]
  query_embeddings <- Embeddings(integrated, reduction = emb_key)[which(integrated[[batch_key]] == query), ]
  nn <- RANN::nn2(data = ref_embeddings, query = query_embeddings, k = k)

  nn_anno <- matrix(anno_tenx[as.vector(nn$nn.idx)], nrow = nrow(nn$nn.idx), ncol = ncol(nn$nn.idx))

  nn_table <- apply(nn_anno, 1, table)

  pred_anno_SS2 <- as.character(lapply(nn_table, function(x) names(x)[which.max(x)]))

  names(pred_anno_SS2) <- colnames(integrated)[which(integrated[[batch_key]] == query)]

  return(pred_anno = pred_anno_SS2)
}

#' Match score between query and reference neighbourhoods
#'
#' Compute, for each query cell, the ratio of the average within-query distance
#' to the average query-to-reference distance in the integrated PCA space.
#' Smaller values imply better alignment.
#'
#' @param integrated  \link[Seurat]{Seurat} object containing the joint embedding
#'   and meta-data.
#' @param ref  Name of the reference dataset/batch (default \code{"10X"}).
#' @param query  Name of the query dataset/batch (default \code{"SS2"}).
#' @param label_key  Name of the meta-data column that stores cell-type / label
#'   information (default \code{"annotation"}).
#' @param emb_key  Name of the reduction to use (default \code{"pca"}).
#' @param batch_key  Name of the meta-data column that indicates dataset/batch
#'   origin (default \code{"dataset"}).
#' @param k  Number of nearest neighbours to consider (default \code{10}).
#'
#' @return
#' Named numeric vector of length equal to the number of query cells.
#'
#' @rdname match_score
#' @export
match_score <- function(integrated, ref = "10X", query = "SS2", label_key = "annotation",
                        emb_key = "pca", batch_key = "dataset", k = 10){
  set.seed(0)

  anno_tenx <- integrated[[label_key]][which(integrated[[batch_key]] == ref)]

  ref_embeddings <- Embeddings(integrated, reduction = emb_key)[which(integrated[[batch_key]] == ref), ]
  query_embeddings <- Embeddings(integrated, reduction = emb_key)[which(integrated[[batch_key]] == query), ]
  nn <- RANN::nn2(data = ref_embeddings, query = query_embeddings, k = k)

  nn_SS2 <- RANN::nn2(data = query_embeddings, query = query_embeddings, k = k)

  metric <- rowMeans(nn_SS2$nn.dists[, -1]) / rowMeans(nn$nn.dists)

  names(metric) <- rownames(query_embeddings)

  return(metric = metric)
}

#' Label transfer via k-NN (soft assignment / probabilities)
#'
#' Identical to \code{\link{label_trans}} but returns a probability matrix
#' instead of a hard label.
#'
#' @param integrated  \link[Seurat]{Seurat} object containing the joint embedding
#'   and meta-data.
#' @param ref  Name of the reference dataset/batch (default \code{"10X"}).
#' @param query  Name of the query dataset/batch (default \code{"SS2"}).
#' @param label_key  Name of the meta-data column that stores cell-type / label
#'   information (default \code{"annotation"}).
#' @param emb_key  Name of the reduction to use (default \code{"pca"}).
#' @param batch_key  Name of the meta-data column that indicates dataset/batch
#'   origin (default \code{"dataset"}).
#' @param k  Number of nearest neighbours to consider (default \code{10}).
#'
#' @return
#' Numeric matrix (query-cells × unique-labels).  Row sums equal 1.
#' Row names are query cell barcodes, column names are reference labels.
#'
#' @rdname label_trans_prob
#' @export
label_trans_prob <- function(integrated, ref = "10X", query = "SS2", label_key = "annotation",
                             emb_key = "pca", batch_key = "dataset", k = 10){
  set.seed(0)

  anno_tenx <- integrated[[label_key]][which(integrated[[batch_key]] == ref)]

  ref_embeddings <- Embeddings(integrated, reduction = emb_key)[which(integrated[[batch_key]] == ref), ]
  query_embeddings <- Embeddings(integrated, reduction = emb_key)[which(integrated[[batch_key]] == query), ]

  nn <- RANN::nn2(data = ref_embeddings, query = query_embeddings, k = k)

  nn_anno <- matrix(anno_tenx[as.vector(nn$nn.idx)], nrow = nrow(nn$nn.idx), ncol = ncol(nn$nn.idx))

  nn_table <- apply(nn_anno, 1, table)

  nn_prob <- lapply(nn_table, function(x) x/sum(x))

  anno_SS2_prob <- matrix(0, nrow = length(nn_prob), ncol = length(unique(anno_tenx)))

  rownames(anno_SS2_prob) <- rownames(query_embeddings)
  colnames(anno_SS2_prob) <- unique(anno_tenx)

  for (i in 1:length(nn_prob)){
    anno_SS2_prob[i, names(nn_prob[[i]])] <- nn_prob[[i]]
  }

  return(anno_SS2_prob = anno_SS2_prob)
}

