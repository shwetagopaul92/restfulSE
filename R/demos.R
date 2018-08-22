#' A set of mouse cortex marker genes.
#' @note \url{http://www.nature.com/doifinder/10.1038/nn.4216}, Fig 1C
#' @examples
#' head(tasicCortex())
#' @return data.frame with columns SYMBOL, GENEID
#' @export
tasicCortex = function() structure(list(SYMBOL = c("Snap25", "Gad1", "Vip", "Sst", "Pvalb", 
                                                   "Slc17a7", "Rorb", "Foxp2", "Aqp4", "Pdgfra", "Mog", "Itgam", 
                                                   "Bgn"), GENEID = c("ENSMUSG00000027273", "ENSMUSG00000070880", 
                                                                      "ENSMUSG00000019772", "ENSMUSG00000004366", "ENSMUSG00000005716", 
                                                                      "ENSMUSG00000070570", "ENSMUSG00000036192", "ENSMUSG00000029563", 
                                                                      "ENSMUSG00000024411", "ENSMUSG00000029231", "ENSMUSG00000076439", 
                                                                      "ENSMUSG00000030786", "ENSMUSG00000031375")), .Names = c("SYMBOL", 
                                                                                                                               "GENEID"), row.names = c(NA, -13L), class = "data.frame")
#' Convenience functions using EC2 server to extract tenx neurons full or subset data
#' @rawNamespace import(rhdf5client, except = groups)
#' @importFrom utils data
#' @importFrom AnnotationDbi keys
#' @importFrom rhdf5client H5S_Array
#' @importFrom S4Vectors SimpleList
#' @param url server URL
#' @param tag string giving the internal dataset name
#' @return RESTfulSummarizedExperiment
#' @examples
#' ss = se100k()
#' # get a set of genes from Tasic et al. 2016 Nature Neuroscience
#' tc = tasicCortex()
#' adultCort = tc$GENEID
#' # subset
#' csums = apply(assay(ss[adultCort,1:500]),1,sum)
#' names(csums) = tc$SYMBOL
#' csums
#' @export
se100k = function(url="http://h5s.channingremotedata.org:5000",
   tag="tenx_100k_sorted") {
  ds = H5S_Array(url, tag)
  ehub = ExperimentHub::ExperimentHub()
  tag = names(AnnotationHub::query(ehub, "st100k"))
  st100k = ehub[[tag[1]]]
  assays(st100k) = SimpleList(counts=ds)
  st100k
}
#' @rdname se100k
#' @aliases se1.3M
#' @note se1.3M provides access to the full 1.3 million neurons
#' with features in their order as given in the original HDF5
#' while se100k provides access to only 100k neurons with
#' expression features sorted by genomic location
#' @return SummarizedExperiment instance
#' @export
se1.3M = function(url="http://h5s.channingremotedata.org:5000",
   tag="tenx_full") {
  ds = H5S_Array(url, tag)
  ehub = ExperimentHub::ExperimentHub()
  tag = names(AnnotationHub::query(ehub, "full_1Mneurons"))
  full_1Mneurons = ehub[[tag[1]]]
  assays(full_1Mneurons) = SimpleList(counts=ds)
  full_1Mneurons
}

#' Convenience function for access to gene-level GTEx tissues, as quantified in recount
#' @import ExperimentHub
#' @param url ip address/host for HDF5 server
#' @param tag name of hdf5 file on server
#' @return SummarizedExperiment instance
#' @examples
#' gtexTiss()
#' @export
gtexTiss = function(url="http://h5s.channingremotedata.org:5000",
   tag="tissues") {
  ds = H5S_Array(url, tag)
  ehub = ExperimentHub::ExperimentHub()
  tag = names(AnnotationHub::query(ehub, "gtexRecount"))
  gtexTiss = ehub[[tag[1]]]
  assays(gtexTiss) = SimpleList(recount=ds)
  gtexTiss
}

#library(restfulSE)
#tiss = gtexTiss()
#binds = grep("Brain", tiss$smtsd); table(tiss$smtsd[binds])

#' Create a data.frame with ENSEMBL and SYMBOL identifiers associated with a GO TERM specified by a regular expression in \code{termPattern}
#' @import GO.db
#' @param termPattern a character string encoding a regular expression to be matched to keys of type TERM in GO.db
#' @param targets \code{columns} to be returned from org.[organism].[inst].db
#' @param organism two-letter code for organism in the OrgDb family of packages
#' @param inst two- or three-letter code (e.g., \code{eg} for ENTREZ GENE or \code{sgd} for yeastgenome.org) identifying institute responsible for annotation
#' @return data.frame
#' @examples
#' gp = goPatt()
#' dim(gp)
#' head(gp)
#' @export
goPatt = function(termPattern="neurotro", 
   targets=c("ENSEMBL", "SYMBOL"), organism="Hs", inst="eg") {
requireNamespace("GO.db")
require(opackn <- paste0("org.", organism, ".", inst, ".db"), character.only=TRUE)
tms = keys(GO.db, keytype="TERM")
nterms = grep(termPattern, tms, value=TRUE) 
nids = AnnotationDbi::select(GO.db, keys=nterms, keytype="TERM", columns=c("GOID", "TERM")) 
AnnotationDbi::select(get(opackn), keys=nids[[2]], 
   keytype="GO", columns=c("ENSEMBL", "SYMBOL"))
}
