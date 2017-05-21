#' convenience function using 2017 Channing server to extract 100k tenx neurons subset
#' @param url server URL
#' @param tag string giving the internal dataset name
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
