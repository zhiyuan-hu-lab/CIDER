#' Compute IDER-Based Similarity
#'
#' Calculate the similarity matrix based on Inter-group Differential Expression 
#' (IDER) metrics with the selected batch effects regressed out.
#'
#' @param seu A Seurat S4 object that includes an \code{initial_cluster} column in its \code{meta.data}. Required.
#' @param group.by.var Character string specifying the column in \code{seu@meta.data} 
#' that defines initial clusters (batch-specific groups). Default is "initial_cluster".
#' @param batch.by.var Character string specifying the metadata column that indicates batch information. Default is "Batch".
#' @param verbose Logical. If \code{TRUE}, progress messages and a progress bar are displayed. Default is \code{TRUE}.
#' @param use.parallel Logical. If \code{TRUE}, parallel computation is used
#'  (requires \code{doParallel}); in this case, no progress bar will be shown. Default is \code{FALSE}.
#' @param n.cores Numeric. The number of cores to use for parallel computing. Default is 1.
#' @param downsampling.size Numeric. The number of cells representing each group. Default is 40.
#' @param downsampling.include Logical. Whether to include groups with fewer 
#' cells than \code{downsampling.size}. Default is \code{FALSE}.
#' @param downsampling.replace Logical. Whether to sample with replacement if a 
#' group is smaller than \code{downsampling.size}. Default is \code{FALSE}.
#'
#' @return A list of objects: a similarity matrix, a numeric vector recording 
#' the cells used, and a data frame of the group combinations included.
#'
#' @seealso \code{\link{plotNetwork}}, \code{\link{finalClustering}}
#'
#' @export
#'
#' @import limma edgeR foreach utils doParallel
#' @importFrom parallel detectCores stopCluster makeCluster
#' @importFrom stats model.matrix cor
#' @importFrom edgeR cpm
getIDEr <- function(seu,
                    group.by.var = "initial_cluster",
                    batch.by.var = "Batch",
                    verbose = TRUE,
                    use.parallel = FALSE,
                    n.cores = 1,
                    downsampling.size = 40,
                    downsampling.include = TRUE,
                    downsampling.replace = TRUE) {

  if(!group.by.var %in% colnames(seu@meta.data)) {
    warning("group.by.var is not in the colnames of seu@meta.data.")
    return(NULL)
  }

  if(!batch.by.var %in% colnames(seu@meta.data)) {
    warning("batch.by.var is not in the colnames of seu@meta.data.")
    return(NULL)
  }

  tmp <- seu@meta.data[,colnames(seu@meta.data) == group.by.var]

  if(batch.by.var != "Batch") { # which is the batch var
    seu$Batch <- seu@meta.data[,colnames(seu@meta.data) == batch.by.var]
  }

  ## merge seurat list
  metadata <- data.frame(
    label = tmp,
    batch = seu$Batch, # batch
    stringsAsFactors = FALSE
  )

  select <- downsampling( # sampling
    metadata = metadata,
    n.size = downsampling.size,
    include = downsampling.include,
    replace = downsampling.replace
  )

  matrix <- as.matrix(.getCountsMatrix(seu)[, select])
  # matrix for dist calculation
  colnames(matrix) <- paste0(colnames(matrix), seq_len(ncol(matrix)))
  # avoid duplication
  keep <- rowSums(matrix > 0.5) > 5
  dge <- edgeR::DGEList(counts = matrix[keep, , drop = FALSE])
  # make a edgeR object
  dge <- dge[!grepl("ERCC-", rownames(dge)), ] # remove ERCC
  dge <- dge[!grepl("MT-", rownames(dge)), ] # remove mitochondria genes
  dge <- dge[!grepl("mt-", rownames(dge)), ]
  # remove mitochondria genes for other species
  logCPM <- cpm(dge, log = TRUE, prior.count = 3) # calculate cpm

  df <- data.frame(
    g = metadata$label[select], ## label
    b = metadata$batch[select], ## batch
    stringsAsFactors = FALSE
  )

  df$detrate <- scale(colMeans(matrix > 0))[, 1] # gene detection rate
  rownames(df) <- colnames(matrix)
  rm(matrix)
  gc()

  N <- length(unique(df$g)) # number of groups

  # get the dataframe of combinations/pairs for comparison
  combinations <- data.frame(g1 = rep(unique(df$g), each = N),
                             g2 = rep(unique(df$g), N),
                             stringsAsFactors = FALSE)
  combinations <- combinations[combinations$g1 != combinations$g2, ]
  combinations$b1 <- df$b[match(combinations$g1, df$g)]
  combinations$b2 <- df$b[match(combinations$g2, df$g)]
  combinations <- combinations[combinations$b1 != combinations$b2, ]

  idx <- c()
  for (i in 2:nrow(combinations)) {
    if (!combinations$g2[i] %in% combinations$g1[seq_len(i - 1)]) {
      idx <- c(idx, i)
    }
  }
  combinations <- combinations[c(1, idx), ]
  rownames(combinations) <- seq_len(nrow(combinations))

  dist_coef <- matrix(0, nrow = N, ncol = N) # distance matrix
  colnames(dist_coef) <- rownames(dist_coef) <- unique(df$g)

  if(use.parallel & n.cores == 1){
    n.cores <- detectCores(logical = FALSE) - 1
    if(verbose){
      message(sprintf("Use %d cores for calculation.", n.cores))
    }
  }

  if (use.parallel == FALSE | n.cores == 1) { # not using parallel -----

    # create progress bar
    if (verbose) {
      message("Generating distance matrix...")
      pb <- txtProgressBar(min = 0, max = nrow(combinations), style = 3)
      k <- 1
    }

    for (i in seq_len(nrow(combinations))) {
      if (verbose) {
        setTxtProgressBar(pb, k) # progress bar
        k <- k + 1
      }

      df$tmp <- "bg"
      df$tmp[which(df$g == combinations$g1[i])] <- "g1"
      df$tmp[which(df$g == combinations$g2[i])] <- "g2"

      design <- model.matrix(~ 0 + tmp + b + detrate, data = df)
      contrast_m <- limma::makeContrasts(
        contrasts = c("tmpg1-tmpbg", "tmpg2-tmpbg"),
        levels = design
      )
      idx1 <- rownames(dist_coef) == combinations$g1[i]
      idx2 <- colnames(dist_coef) == combinations$g2[i]
      dist_coef[idx1, idx2] <- getGroupFit(logCPM, design, contrast_m)

    }
    if (verbose == TRUE) {
      close(pb) # close progress bar
    }
  } else if (use.parallel == TRUE) { # using parallel -----

    # decide OS and register parallel the parallel backend
    if(Sys.info()["sysname"] == "Linux") {
      registerDoParallel(cores = n.cores)
    } else{
      if(n.cores > detectCores(logical = FALSE) - 1){
        warning("too many cores assign. Setting n.cores =
                detectCores(logical = FALSE) - 1")
        n.cores <- detectCores(logical = FALSE) - 1
      }
      cl <- makeCluster(n.cores)
      registerDoParallel(cl = cl)
    }

    n.iter <- nrow(combinations) # number of iterations

    i <- NULL
    j <- NULL


    df_dist <- foreach(i = combinations$g1, j = combinations$g2,
                       df = rep(list(df), n.iter),
                       logCPM = rep(list(logCPM), n.iter),
                       .combine = "rbind", .verbose = verbose) %dopar%
                       {

                         df$tmp <- "bg"
                         df$tmp[df$g == i] <- "g1"
                         df$tmp[df$g == j] <- "g2"

                         design <- model.matrix(~  0 + tmp + b + detrate,
                                                data = df)
                         contrast_m <- limma::makeContrasts(
                           contrasts = c("tmpg1-tmpbg", "tmpg2-tmpbg"),
                           levels = design
                         )
                         getGroupFit(logCPM, design, contrast_m)
                       }

    if(Sys.info()["sysname"] == "Linux") {
      stopImplicitCluster()
    } else {
      stopCluster(cl)
    }

    for (i in seq_len(nrow(combinations))) {
      # assign the distance values into the matrix
      idx1 <- rownames(dist_coef) == combinations$g1[i]
      idx2 <- colnames(dist_coef) == combinations$g2[i]
      dist_coef[idx1, idx2] <- df_dist[i, 1]
    }
  }

  return(list(dist_coef + t(dist_coef), select, combinations))
}

#' Calculate IDER-Based Similarity Between Two Groups
#'
#' This function calculates the IDER-based similarity between two groups using a linear model.
#'
#' @param logCPM A numeric matrix of log-transformed counts per million.
#' @param design A design matrix for the differential expression analysis.
#' @param contrast_m A contrast matrix specifying the comparison between the two groups.
#'
#' @return A numeric value representing the IDER-based similarity between the two groups.
#'
#' @importFrom stats cor .lm.fit
getGroupFit <- function(logCPM, design, contrast_m){
  fit <- .lm.fit(design, t(logCPM))
  coef <- .zeroDominantMatrixMult(t(fit$coefficients), contrast_m)
  return(cor(coef[, 1], coef[, 2]))
}

