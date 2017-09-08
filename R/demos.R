#' convenience functions using EC2 server to extract tenx neurons full or subset data
#' @import rhdf5client
#' @importFrom AnnotationDbi select keys
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
se100k = function(url="http://54.174.163.77:5000",
   tag="tenx_100k_sorted") {
  src = H5S_source(url)
  ds = src[[tag]]
  data(st100k)
  RESTfulSummarizedExperiment(st100k, ds)
}
#' @rdname se100k
#' @aliases se1.3M
#' @note se1.3M provides access to the full 1.3 million neurons
#' with features in their order as given in the original HDF5
#' while se100k provides access to only 100k neurons with
#' expression features sorted by genomic location
#' @export
se1.3M = function(url="http://54.174.163.77:5000",
   tag="tenx_full") {
  src = H5S_source(url)
  ds = src[[tag]]
  data(full_1Mneurons)
  full_1Mneurons = as(full_1Mneurons, "RangedSummarizedExperiment")
  RESTfulSummarizedExperiment(full_1Mneurons, ds)
}

#' convenience function for access to gene-level GTEx tissues, as quantified in recount
#' @param url ip address/host for HDF5 server
#' @param tag name of hdf5 file on server
#' @return RESTfulSummarizedExperiment instance
#' @examples
#' gtexTiss()
#' @export
gtexTiss = function(url="http://54.174.163.77:5000",
   tag="tissues") {
  src = H5S_source(url)
  ds = src[[tag]]
  data(gtexRecount)
  gtexTiss = as(gtexRecount, "RangedSummarizedExperiment")
  RESTfulSummarizedExperiment(gtexTiss, ds)
}

#library(restfulSE)
#tiss = gtexTiss()
#binds = grep("Brain", tiss$smtsd); table(tiss$smtsd[binds])

#' create a data.frame with ENSEMBL and SYMBOL identifiers associated with a GO TERM specified by a regular expression in \code{termPattern}
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
