#' hdf5server-based assay for SummarizedExperiment
#' @import SummarizedExperiment
#' @exportClass RESTfulSummarizedExperiment
setClass("RESTfulSummarizedExperiment",
   contains="RangedSummarizedExperiment", 
     representation(source="H5SDatasets",
                    globalDimnames="list"))

#' construct RESTfulSummarizedExperiment
#' @param se SummarizedExperiment instance, assay component can be empty SimpleList
#' @param source instance of H5SDatasets
#' @examples
#' \dontrun{
#' bb = banoH5()
#' data(banoSEMeta)
#' rr = RESTfulSummarizedExperiment(banoSEMeta, bb)
#' rr
#' assay(rr)
#' }
#' @export RESTfulSummarizedExperiment
RESTfulSummarizedExperiment = function(se, source) {
   stopifnot(is(se, "RangedSummarizedExperiment")) # for now
   new("RESTfulSummarizedExperiment", se, source=source,
        globalDimnames=dimnames(se))
}

setMethod("assayNames", "RESTfulSummarizedExperiment", function(x, ...) {
 "(served by HDF5Server)"
})

#' @exportMethod [
setMethod("[", c("RESTfulSummarizedExperiment",
     "numeric", "numeric", "ANY"), function(x,i,j,...,drop=FALSE) {
   x = BiocGenerics:::replaceSlots(x, rowRanges = rowRanges(x)[i],
                         colData = colData(x)[j,],
                         check=FALSE)
   x
   })

setMethod("assays", c("RESTfulSummarizedExperiment"), function(x, ...,
   withDimnames=TRUE) {
#   warning("use assay(), only one allowed at present for RESTful SE")
#   assay(x, ...)  # document properly
   SimpleList("placeholder")
})
 
#' @exportMethod assay
setMethod("assay", c("RESTfulSummarizedExperiment", 
      "missing"), function(x, i, ...) {
       rowsToGet = match(rownames(x), x@globalDimnames[[1]])
       colsToGet = match(colnames(x), x@globalDimnames[[2]])
       ans = x@source[rowsToGet, colsToGet]
       dimnames(ans) = dimnames(x)
       ans
})
