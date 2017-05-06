#' banoSEMeta: metadata RangedSummarizedExperiment shell for banovichSE
#' @docType data
#' @format RangedSummarizedExperiment instance
#' @source 10.1371/journal.pgen.1004663
"banoSEMeta"
#' gr450k: GRanges with metadata for illumina 450k methylation assay
#' @docType data
#' @format RangedSummarizedExperiment instance
#' @source Bioconductor
"gr450k"
#' tenx_100k_sorted: metadata RangedSummarizedExperiment shell for 100k cells from 10x genomics mouse 1.3 million cell dataset
#' @docType data
#' @format RangedSummarizedExperiment instance
#' @source \url{https://community.10xgenomics.com/t5/10x-Blog/Our-1-3-million-single-cell-dataset-is-ready-to-download/ba-p/276}
#' @note used github.com/mtmorgan/TENxGenomics package to obtain 
#' SummarizedExperiment, then added range information, sorted within
#' chromosome, and saved shell for use with HDF5 server
"tenx_100k_sorted"
