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
#' st100k: metadata RangedSummarizedExperiment shell for 100k cells from 10x genomics 1.3 million neuron dataset
#' @docType data
#' @format RangedSummarizedExperiment instance
#' @source \url{https://community.10xgenomics.com/t5/10x-Blog/Our-1-3-million-single-cell-dataset-is-ready-to-download/ba-p/276}
#' @note used github.com/mtmorgan/TENxGenomics package to obtain 
#' SummarizedExperiment, then added range information, sorted within
#' chromosome, and saved shell for use with HDF5 server
"st100k"

#' A set of mouse cortex marker genes.
#' @note \url{http://www.nature.com/doifinder/10.1038/nn.4216}, Fig 1C
#' @export
tasicCortex = function() structure(list(SYMBOL = c("Snap25", "Gad1", "Vip", "Sst", "Pvalb", 
"Slc17a7", "Rorb", "Foxp2", "Aqp4", "Pdgfra", "Mog", "Itgam", 
"Bgn"), GENEID = c("ENSMUSG00000027273", "ENSMUSG00000070880", 
"ENSMUSG00000019772", "ENSMUSG00000004366", "ENSMUSG00000005716", 
"ENSMUSG00000070570", "ENSMUSG00000036192", "ENSMUSG00000029563", 
"ENSMUSG00000024411", "ENSMUSG00000029231", "ENSMUSG00000076439", 
"ENSMUSG00000030786", "ENSMUSG00000031375")), .Names = c("SYMBOL", 
"GENEID"), row.names = c(NA, -13L), class = "data.frame")

