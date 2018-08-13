# isplit/sproc from rhdf5client

#' HDF5Server-based assay for SummarizedExperiment
#' @import SummarizedExperiment
#' @importClassesFrom rhdf5client H5S_dataset
#' @exportClass RESTfulSummarizedExperiment
setClass("RESTfulSummarizedExperiment",
   contains="RangedSummarizedExperiment", 
     representation(source="H5S_dataset",
                    globalDimnames="list"))

#' Construct RESTfulSummarizedExperiment
#' @aliases RESTfulSummarizedExperiment,RangedSummarizedExperiment,H5S_dataset-method
#' @param se SummarizedExperiment instance, assay component can be empty SimpleList
#' @param source instance of H5S_dataset
#' @examples
#' require("rhdf5client")
#' bigec2 = H5S_source(serverURL="http://h5s.channingremotedata.org:5000")
#' banoh5 = bigec2[["assays"]] # banovichSE
#' ehub = ExperimentHub::ExperimentHub()
#' myfiles <- AnnotationHub::query(ehub , "restfulSEData")
#' myfiles[["EH551"]] -> banoSEMeta
#' rr = RESTfulSummarizedExperiment(banoSEMeta, banoh5)
#' rr
#' rr2 = rr[1:4, 1:5] # just modify metadata
#' rr2
#' assay(rr2) # extract data
#' @return instance of RESTfulSummarizedExperiment
#' @exportMethod RESTfulSummarizedExperiment
#' @export RESTfulSummarizedExperiment
setGeneric("RESTfulSummarizedExperiment",
  function(se, source) standardGeneric("RESTfulSummarizedExperiment"))
setMethod("RESTfulSummarizedExperiment", c("RangedSummarizedExperiment",
   "H5S_dataset"), function(se, source) {
 .RESTfulSummarizedExperiment(se, source)
})

#' hidden constructor
#' @rdname RESTfulSummarizedExperiment
#' @rawNamespace importFrom("methods", "as", "is", "new", "slot", "slot<-", "validObject")
#' @importFrom rhdf5client internalDim
.RESTfulSummarizedExperiment = function(se, source) {
   stopifnot(is(se, "RangedSummarizedExperiment")) # for now
   d = internalDim(source)
   if (!all(d == rev(dim(se)))) {
       cat("rev(internal dimensions of H5S_dataset) is", rev(d), "\n")
       cat("dim(se) is", dim(se), "\n")
       stop("these must agree.\n")
       }
   new("RESTfulSummarizedExperiment", se, source=source,
        globalDimnames=dimnames(se))
}

#' @rdname RESTfulSummarizedExperiment-class
#' @aliases assayNames,RESTfulSummarizedExperiment-method
setMethod("assayNames", "RESTfulSummarizedExperiment", function(x, ...) {
 "(served by HDF5Server)"
})

#' @rdname RESTfulSummarizedExperiment-class
#' @importFrom DelayedArray rowRanges
#' @aliases [,RESTfulSummarizedExperiment,numeric,numeric,ANY-method
#' @param x instance of RESTfulSummarizedExperiment
#' @param i numeric selection vector
#' @param j numeric selection vector
#' @param \dots not used
#' @param drop not used
#' @return instance of RESTfulSummarizedExperiment
#' @exportMethod [
setMethod("[", c("RESTfulSummarizedExperiment",
     "numeric", "numeric", "ANY"), function(x,i,j,...,drop=FALSE) {
  if (is(x, "RangedSummarizedExperiment")) {
   x = replaceSlots(x, rowRanges = rowRanges(x)[i],
                         colData = colData(x)[j,],
                         check=FALSE)
   }
  else if (is(x, "SummarizedExperiment")) {
   x = replaceSlots(x, rowData = rowData(x)[i],
                         colData = colData(x)[j,],
                         check=FALSE)
   }
   x
   })

#' @name assay
#' @rdname RESTfulSummarizedExperiment
#' @importFrom rhdf5client isplit sproc
#' @note RESTfulSummarizedExperiment contains a global dimnames
#' list generated at creation.  It is possible that standard operations 
#' on a SummarizedExperiment will engender dimnames components that
#' differ from the initial global dimnames, principally through
#' uniqification (adding suffixes when dimname elements are
#' repeated).  When this is detected, assay() will fail with a complaint
#' about length(setdiff(*names(x), x@globalDimnames[[...]])).
#' @aliases assay,RESTfulSummarizedExperiment,missing-method
#' @param x instance of RESTfulSummarizedExperiment
#' @param i not used
#' @param \dots not used
#' @return matrix
#' @exportMethod assay
setMethod("assay", c("RESTfulSummarizedExperiment", "missing"), 
    function(x, i, ...) {
    stopifnot(length(rownames(x))>0)
    stopifnot(length(colnames(x))>0)
    stopifnot(length(setdiff(rownames(x), x@globalDimnames[[1]]))==0)
    stopifnot(length(setdiff(colnames(x), x@globalDimnames[[2]]))==0)
    rowsToGet = match(rownames(x), x@globalDimnames[[1]])
    colsToGet = match(colnames(x), x@globalDimnames[[2]])
    ind1 = sproc(isplit(colsToGet))  # may need to be double loop
    ind2 = sproc(isplit(rowsToGet))
#    if (length(ind1)>1 | length(ind2)>1) warning("as of 5/5/17 only processing contiguous block requests, will generalize soon; using first block only")
    if (length(ind1)==1 & length(ind2)==1) 
       ans = t(x@source[ ind1[[1]], ind2[[1]] ])
    else if (length(ind2)==1) {
       ansl = lapply(ind1, function(i1) t(x@source[i1, ind2[[1]] ]))
       ans = do.call(cbind,ansl)
       }
    else if (length(ind1)==1) {
       ansl = lapply(ind2, function(i2) t(x@source[ind1[[1]], i2 ]))
       ans = do.call(rbind,ansl)
       }
    else {
       ansl = lapply(ind1, function(i1) 
                do.call(rbind, lapply(ind2, 
                  function(i2) t(x@source[i1, i2]))))
       ans = do.call(cbind, ansl)
         }


    dimnames(ans) = list(x@globalDimnames[[1]][rowsToGet], 
                x@globalDimnames[[2]][colsToGet])
    ans
})

#' Assays access for RESTfulSummarizedExperiment
#' @importFrom S4Vectors SimpleList
#' @param x instance of RESTfulSummarizedExperiment
#' @param \dots not used
#' @param withDimnames logical defaults to TRUE
#' @return at present a SimpleList is returned as a dummy placeholder
#' @exportMethod assays
setMethod("assays", c("RESTfulSummarizedExperiment"), function(x, ...,
   withDimnames=TRUE) {
#   warning("use assay(), only one allowed at present for RESTful SE")
#   assay(x, ...)  # document properly
   SimpleList("placeholder")
})
 

#' Dimension access for RESTfulSummarizedExperiment
#' @param x instance of RESTfulSummarizedExperiment
#' @return vector of nrows, ncols
#' @exportMethod dim
setMethod("dim", "RESTfulSummarizedExperiment", function(x)
   c(length(rownames(x)), length(colnames(x)))
)
