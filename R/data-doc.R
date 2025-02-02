#' Dendritic Data
#'
#' A Seurat object containing dendritic cell data.
#'
#' @description This dataset comprises data for four dendritic cell subtypes (CD141, 
#' CD1C, DoubleNeg, and pDC) from two batches. The raw count matrix and sample 
#' information were downloaded from a curated set, and cells with fewer than 
#' 500 detected genes have been removed.
#'
#' @format A Seurat object.
#'
#' @usage data(dendritic)
#'
#' @source The data were downloaded from \url{https://hub.docker.com/r/jinmiaochenlab/batch-effect-removal-benchmarking}.  
#' Reference 1: Tran HTN, Ang KS, Chevrier M, Lee NYS, Goh M, Chen J. A benchmark of batch-effect correction methods for single-cell RNA sequencing data. Genome Biol. 2020;21:1–32.  
#' Reference 2: Villani A-C, Satija R, Reynolds G, Sarkizova S, Shekhar K, Fletcher J, et al. Single-cell RNA-seq reveals new types of human blood dendritic cells, monocytes, and progenitors. Science. 2017;356(6335).
#'
"dendritic"

#' Pancreas Count Matrix
#'
#' A sparse count matrix for cross-species pancreatic data.
#'
#' @description This dataset comprises a single-cell RNA-seq count matrix for pancreatic data from human and mouse samples.
#'
#' @format A sparse matrix of class \code{"dgCMatrix"} with 12474 rows (genes) and 10127 columns (cells).
#'
#' @usage data(pancreas_counts)
#'
#' @source The count matrix and sample information were downloaded from NCBI GEO accession GSE84133.  
#' Reference: Baron M, Veres A, Wolock SL, Faust AL, Gaujoux R, Vetere 
#' A, et al. A single-cell transcriptomic map of the human and mouse pancreas reveals inter- 
#' and intra-cell population structure. Cell Syst. 2016;3:346–360.e4.
#'
"pancreas_counts"

#' Pancreas Metadata
#'
#' Cell-level metadata for cross-species pancreatic data.
#'
#' @description This dataset provides cell-level metadata for the human and mouse 
#' pancreatic data used in the study.
#'
#' @format A data frame with 10127 rows and 3 columns:
#' \describe{
#'   \item{Batch}{Species information (human or mouse).}
#'   \item{Group}{Cell type annotation.}
#'   \item{Sample}{Donor information.}
#' }
#'
#' @usage data(pancreas_meta)
#'
#' @source The metadata were downloaded alongside the count matrix from NCBI GEO accession GSE84133.  
#' Reference: Baron M, Veres A, Wolock SL, Faust AL, Gaujoux R, Vetere A, et al. 
#' A single-cell transcriptomic map of the human and mouse pancreas reveals inter- and intra-cell 
#' population structure. Cell Syst. 2016;3:346–360.e4.
#'
"pancreas_meta"

