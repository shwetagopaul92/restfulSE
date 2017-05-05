
# utilities for index processing
# sproc(isplit(vec)) will convert vec representing R integer vector
# into a list of HDF5server 'select' index candidates
#' isplit converts a numeric vector into a list of sequences for compact reexpression
#' @name isplit
#' @rdname sproc
#' @param x a numeric vector (should be integers)
#' @export
isplit = function(x) {
 dx = diff(x)
 rdx = rle(dx)
 grps = c(1, rep(1:length(rdx$length), rdx$length))
 split(x, grps)
}

#' sproc massages output of isplit into HDF5 select candidates
#' @name sproc
#' @rdname sproc
#' @param spl output of isplit
#' @examples
#' inds = c(1:10, seq(25,50,2), seq(200,150,-2))
#' sproc(isplit(inds))
#' @export
sproc = function(spl) {
# spl is output of isplit
ans = lapply(spl, function(x) {
   if (length(x)==1) return(paste(x-1,":",x,":1", sep=""))
   d = x[2]-x[1]
   return(paste(x[1]-1, ":", x[length(x)], ":", as.integer(d),
     sep=""))
   })
ans
}

#myvec = myvec = c(2:6, 12, 17, seq(30,7,-2))
#sproc(isplit(myvec))

#' hdf5server-based assay for SummarizedExperiment
#' @import SummarizedExperiment
#' @exportClass RESTfulSummarizedExperiment
setClass("RESTfulSummarizedExperiment",
   contains="RangedSummarizedExperiment", 
     representation(source="H5S_dataset",
                    globalDimnames="list"))

#' construct RESTfulSummarizedExperiment
#' @param se SummarizedExperiment instance, assay component can be empty SimpleList
#' @param source instance of H5S_dataset
#' @examples
#' bigec2 = H5S_source(serverURL="http://54.174.163.77:5000")
#' n100k = bigec2[["neurons100k"]]
#' data(tenx_100k_sorted)
#' rr = RESTfulSummarizedExperiment(tenx_100k_sorted, n100k)
#' rr
#' rr2 = rr[1:4, 1:5] # just modify metadata
#' rr2
#' assay(rr2) # extract data
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

#' @name assay
#' @rdname RESTfulSummarizedExperiment
#' @exportMethod assay
setMethod("assay", c("RESTfulSummarizedExperiment", "missing"), 
    function(x, i, ...) {
    rowsToGet = match(rownames(x), x@globalDimnames[[1]])
    colsToGet = match(colnames(x), x@globalDimnames[[2]])
    ind1 = sproc(isplit(colsToGet))  # may need to be double loop
    ind2 = sproc(isplit(rowsToGet))
    if (length(ind1)>1 | length(ind2)>1) warning("as of 5/5/17 only processing contiguous block requests, will generalize soon; using first block only")
    ans = t(x@source[ ind1[[1]], ind2[[1]] ])
    dimnames(ans) = list(x@globalDimnames[[1]][rowsToGet], 
                x@globalDimnames[[2]][colsToGet])
    ans
})

#' @exportMethod assays
setMethod("assays", c("RESTfulSummarizedExperiment"), function(x, ...,
   withDimnames=TRUE) {
#   warning("use assay(), only one allowed at present for RESTful SE")
#   assay(x, ...)  # document properly
   SimpleList("placeholder")
})
 

#' @exportMethod dim
setMethod("dim", "RESTfulSummarizedExperiment", function(x)
   c(length(rownames(x)), length(colnames(x)))
)
