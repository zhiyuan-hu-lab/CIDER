#' Dendritic Data
#'
#' @description The data contained four cell subtypes (CD141, CD1C, DoubleNeg, and pDC) from two batches.
#' The raw count matrix and the sample information were also downloaded from the curated set.
#' Cells with less than 500 genes detected were removed. 
#'
#' @format A Seurat Object
#' @usage data(dendritic)
#'
#' @source The data were downloaded from https://hub.docker.com/r/jinmiaochenlab/batch-effect-removal-benchmarking
#' ref 1: Tran HTN, Ang KS, Chevrier M, Lee NYS, Goh M, Chen J. A benchmark of batch-effect correction methods for single-cell RNA sequencing data. Genome Biol. 2020;21:1–32.
#' ref 2: Villani A-C, Satija R, Reynolds G, Sarkizova S, Shekhar K, Fletcher J, et al. Single-cell RNA-seq reveals new types of human blood dendritic cells, monocytes, and progenitors. Science. 2017;356(6335).
"dendritic"


#' Pancreas Count Matrix
#'
#' @description Single-cell RNA-seq count matrix for human and mouse cross-species pancreatic data
#'
#' @format sparse Matrix of class "dgCMatrix" with 12474 rows (genes) and 10127 columns (cells).
#'
#' @usage data(pancreas_counts)
#'
#' @source The count matrix and sample information were downloaded from NCBI GEO accession GSE84133. 
#' ref: Baron M, Veres A, Wolock SL, Faust AL, Gaujoux R, Vetere A, et al. A single-cell transcriptomic map of the human and mouse pancreas reveals inter- and intra-cell population structure. Cell Syst. 2016;3:346–360.e4.
#'
"pancreas_counts"

#' Pancreas Metadata
#'
#' @description Cell-level metadata for human and mouse cross-species pancreatic data
#'
#' @format A data frame with 10127 rows and 3 columns:
#' \describe{
#'   \item{Batch}{Species information, i.e., human or mouse}
#'   \item{Group}{Cell type annotation}
#'   \item{Sample}{Donor information}
#' }
#'
#' @usage data(pancreas_meta)
#'
#' @source The count matrix and sample information were downloaded from NCBI GEO accession GSE84133. 
#' ref: Baron M, Veres A, Wolock SL, Faust AL, Gaujoux R, Vetere A, et al. A single-cell transcriptomic map of the human and mouse pancreas reveals inter- and intra-cell population structure. Cell Syst. 2016;3:346–360.e4.
#'
"pancreas_meta"
