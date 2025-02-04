#' Initial Clustering for Evaluating Integration
#'
#' This function applies HDBSCAN, a density-based clustering algorithm, to the corrected dimension 
#' reduction of a Seurat object.
#'
#' @param seu A Seurat object containing integrated or batch-corrected data (e.g. PCA results).
#' @param batch.var Character string specifying the metadata column that contains batch information. 
#' Default is "Batch".
#' @param reduction Character string specifying the name of the dimension reduction to use (e.g. "PCA"). 
#' Default is "PCA".
#' @param dims Numeric vector indicating the dimensions to be used for initial clustering. Default is 1:15.
#' @param minPts Integer specifying the minimum number of points required to form a cluster. 
#' This value is passed to the \code{hdbscan} function. Default is 25.
#'
#' @return A Seurat object with two additional columns in its \code{meta.data}: 
#' \code{dbscan_cluster} and \code{initial_cluster}.
#'
#' @seealso \code{\link{getIDEr}}, \code{\link{estimateProb}}
#'
#' @export
#'
#' @import Seurat
#' @importFrom dbscan hdbscan
hdbscan.seurat <- function(seu, batch.var = "Batch", reduction = "pca",
                            dims = seq_len(15), minPts = 25){
  if(!reduction %in% Reductions(seu)) stop("Reduction does not exist. Please check the reduction paramter.")
  seu <- RunTSNE(seu, reduction = reduction, dims = dims, reduction.name = "tsne")
  
  tsne_embeddings <- Embeddings(seu, "tsne")
  # tsne_embeddings <- if (.checkSeuratObjectVersion(seu) == "v5") {
  #   Embeddings(seu, "tsne")
  # } else {
  #   Reductions(seu, "tsne")@cell.embeddings
  # }

  res <- hdbscan(tsne_embeddings[,seq_len(2)],
                  minPts = minPts)
  seu$dbscan_cluster <- factor(as.character(res$cluster))
  seu$initial_cluster <- factor(paste0(seu$dbscan_cluster, "_", seu@meta.data[[batch.var]]))
  return(seu)
}

#' Estimate the Empirical Probability of Whether Two Set of Cells
#' from Distinct Batches Belong to the Same Population
#'
#' This function computes the empirical probability that two sets of cells from 
#' distinct batches belong to the same population, based on the output of \code{getIDEr}.
#'
#' @param seu A Seurat object.
#' @param ider A list returned by the \code{getIDEr} function.
#' @param batch.var Character string specifying the metadata column that contains 
#' batch information. Default is "Batch".
#' @param n_size Numeric value indicating the number of cells per group used to 
#' compute the similarity. Default is 40.
#' @param n.perm Numeric value specifying the number of permutations to perform.
#' @param verbose Logical. If \code{TRUE}, progress messages are printed. 
#' Default is \code{FALSE}.
#'
#' @return A Seurat object with additional columns for the IDER-based similarity 
#' and the empirical probability of rejection.
#'
#' @seealso \code{\link{hdbscan.seurat}}, \code{\link{getIDEr}}
#'
#' @export
#'
#' @import limma edgeR foreach utils doParallel
#' @importFrom kernlab specc
estimateProb <- function(seu, ider, batch.var = "Batch", n_size = 40,
                          #seeds = c(12345, 89465, 10385, 10385, 17396),
                          n.perm = 5, verbose = FALSE){

  dist_coef <- ider[[1]]
  dist_coef[upper.tri(dist_coef)] <- 0
  # Select positive control
  pos_control <- c(rownames(dist_coef)[which.max(apply(dist_coef, 1, max))],
                   colnames(dist_coef)[which.max(apply(dist_coef, 2, max))])

  idx <- seu$initial_cluster %in% pos_control

  combinations_all <- c()
  bg_dist_coef_list <- list() # background distance distribution
  seu_selected <- seu[,idx]
  # pca <- seu@reductions$pca@cell.embeddings[idx, seq_len(15)]
  pca <- Embeddings(seu, "pca")[idx, seq_len(15)]
  # use first 15 PCs for spectral clustering

  for(itor in seq_len(n.perm)) {
    # force the positive control group into first groups
    # set.seed(seeds[itor])
    res <- specc(pca, centers = 5) # spectral clustering

    ## Calculate background distribution
    seu_selected$forced_cluster <- res@.Data
    seu@meta.data$forced_cluster <-
      seu_selected$forced_cluster[match(colnames(seu),
                                        colnames(seu_selected))]
    seu$forced_cluster[is.na(seu$forced_cluster)] <- "bg"
    metadata <- data.frame(label = paste0(seu$forced_cluster, "_", seu@meta.data[[batch.var]]),
                            batch = seu@meta.data[[batch.var]],
                            stringsAsFactors = FALSE)

    # downsampling
    select <- downsampling(metadata = metadata, n.size = n_size,
                           include = TRUE, replace = FALSE, lower.cutoff = 15)

    # IDER
    # matrix <- as.matrix(seu@assays$RNA@counts[,select])
    # version-aware layer/slot handling
    matrix <- .getCountsMatrix(seu)[, select]
    
    keep <- rowSums(matrix > 0.5) > 5
    dge <- edgeR::DGEList(counts = matrix[keep, , drop = FALSE])
    # make a edgeR object
    dge <- dge[!grepl("ERCC-", rownames(dge)),] # remove ERCC
    dge <- dge[!grepl("MT-", rownames(dge)),]

    df <- data.frame(g = metadata$label[select],
                     b = metadata$batch[select], ## batch
                     stringsAsFactors = FALSE) ## label
    df$detrate <- scale(colMeans(matrix > 0))[, 1]
    colnames(matrix) <- paste0(colnames(matrix), "_", seq_len(ncol(matrix)))
    rownames(df) <- colnames(matrix)

    GROUPS <- unique(df$g)
    N <- length(GROUPS)

    combinations <- data.frame(g1 = rep(unique(df$g), each = N),
                               g2 = rep(unique(df$g), N),
                               stringsAsFactors = FALSE)
    combinations <- combinations[combinations$g1 != combinations$g2, ]
    combinations$b1 <- df$b[match(combinations$g1, df$g)]
    combinations$b2 <- df$b[match(combinations$g2, df$g)]
    combinations <- combinations[combinations$b1 != combinations$b2, ]
    idx <- c()
    for(i in 2:nrow(combinations)){
      if(!combinations$g2[i] %in% combinations$g1[seq_len(i-1)]) {
        idx <- c(idx, i)
      }
    }

    combinations <- combinations[c(1,idx),]
    rownames(combinations) <- seq_len(nrow(combinations))
    combinations <-
      combinations[!combinations$g1 %in% c("bg_Batch1", "bg_Batch2")  &
                     !combinations$g2 %in% c("bg_Batch1", "bg_Batch2"),]
    combinations$similarity <- NA
    combinations$iteration <- itor

    bg_dist_coef <- matrix(0, nrow = N, ncol = N)
    colnames(bg_dist_coef) <- rownames(bg_dist_coef) <- GROUPS

    # create progress bar
    if (verbose == TRUE) {
      message("Generating distance matrix...")
      pb <- txtProgressBar(min = 0, max = nrow(combinations), style = 3)
      k <- 1
    }
    logCPM_all <- cpm(dge, log = TRUE, prior.count = 3)
    for (i in seq_len(nrow(combinations))){
      if (verbose == TRUE) {
        setTxtProgressBar(pb, k) # progress bar
        k <- k+1
      }
      df$tmp <- NA
      df$tmp[df$g %in% c("bg_Batch1", "bg_Batch2")] <- "bg"
      df$tmp[df$g == combinations$g1[i]] <- "g1"
      df$tmp[df$g == combinations$g2[i]] <- "g2"

      idx <- which(!is.na(df$tmp))
      design <- model.matrix(~  0 + tmp + b + detrate, data = df[idx, ])
      contrast_m <- makeContrasts(contrasts = c("tmpg1-tmpbg", "tmpg2-tmpbg"),
                                  levels = design)
      logCPM <- logCPM_all[,idx]
      fit <- lmFit(logCPM, design)
      group_fit <- contrasts.fit(fit, contrast_m)

      idx1 <- rownames(bg_dist_coef) == combinations$g1[i]
      idx2 <- colnames(bg_dist_coef) == combinations$g2[i]
      combinations$similarity[i] <- bg_dist_coef[idx1, idx2] <-
        cor(coef(group_fit)[,1], coef(group_fit)[,2])
    }

    combinations_all <- rbind(combinations_all, combinations)
    bg_dist_coef_list[[itor]] <- bg_dist_coef

    if(verbose == TRUE) {
      close(pb) # close progress bar
    }
  }

  idx <- getSharedGroups(seu, ider[[1]], batch.var = batch.var)
  shared_g <- idx[[1]]
  idx1 <- idx[[2]]
  idx2 <- idx[[3]]

  p_mat <- getProbability(ider[[1]][idx1, idx2], combinations_all$similarity)

  # assign similiary
  scores <- diag(ider[[1]][idx1, idx2])
  names(scores) <- shared_g
  seu@meta.data$similarity <- scores[match(seu$dbscan_cluster, names(scores))]

  # assign p values
  scores <- diag(p_mat[idx1, idx2])
  names(scores) <- shared_g
  seu@meta.data$pvalue <- scores[match(seu$dbscan_cluster, names(scores))]

  return(seu)
}

getProbability <- function(x, bg_similarity) { # calculate the probability
  try({
    res <- matrix(NA, ncol = ncol(x), nrow = nrow(x))
    rownames(res) <- rownames(x)
    colnames(res) <- colnames(x)
    len <- length(bg_similarity)
    for(i in seq_len(nrow(x))){
      for(j in seq_len(ncol(x))){
        res[i,j] <- sum(bg_similarity > x[i,j]) / len
      }
    }
    return(res)
  })
}
