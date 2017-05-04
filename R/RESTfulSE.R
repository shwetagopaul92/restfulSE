#' hdf5server-based assay for SummarizedExperiment
#' @import SummarizedExperiment
#' @exportClass RESTfulSummarizedExperiment
setClass("RESTfulSummarizedExperiment",
   contains="RangedSummarizedExperiment", 
     representation(source="RESTfulH5",
                    globalDimnames="list"))

#' construct RESTfulSummarizedExperiment
#' @param se SummarizedExperiment instance, assay component can be empty SimpleList
#' @param source instance of H5SDatasets
#' @examples
#' \dontrun{
#' shi = H5S_source(serverURL="http://170.223.248.164:7248")
#' lin1 = links(shi,1) 
#' ds1 = datasetRefs(lin1,1,drop=1:4)
#' n100k = ds1[["neurons100k"]]
#' data(tenx_100k_sorted)
#' rr = RESTfulSummarizedExperiment(tenx_100k_sorted, n100k)
#' rr
#' assay(rr[1:10,1:20])
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

#' @exportMethod assays
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

#' @exportMethod dim
setMethod("dim", "RESTfulSummarizedExperiment", function(x)
   c(length(rownames(x)), length(colnames(x)))
)
