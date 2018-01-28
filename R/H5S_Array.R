# 
# implementation of 
# https://github.com/Bioconductor/DelayedArray/blob/master/vignettes/02-Implementing_a_backend.Rmd
# for HDF server back end, based on rhdf5client
#
# key contract: if is (X, "H5S_source") and host is the
# name of an HDF5 file on the server, in /data, then Z = X[[ host ]]
# returns an H5S_dataset instance; this in turn answers to Z[i,j]
# with a numeric matrix
#
# longrunning jan 2018
# domain http://54.174.163.77:5000
# host tenx_full

#' H5S_Array for HDF Server content
#' @import DelayedArray
setClass("H5S_ArraySeed",
   contains="Array",
   slots = c(
     filepath="character",
     domain="character",
     host="character",
     H5S_dataset="H5S_dataset"
     ))

H5S_ArraySeed = function(filepath, host) {
  requireNamespace("rhdf5client")
  dom = try(rhdf5client::H5S_source(filepath))
  if (!is(dom, "H5S_source")) stop("could not resolve H5S_source request on filepath")
#
# for HDF Server back end, the following invocation of [[
# establishes the link to numerical data source
#
  ds = try(dom[[host]])
  if (!is(ds, "H5S_dataset")) stop("could not resolve H5S_dataset request on filepath[[host]]")
  new("H5S_ArraySeed", filepath=filepath,
         domain=filepath, host=host, H5S_dataset=ds)
  }

#' dimnames not stored with H5S_source as of Jan 2018
#' @param x instance of H5S_ArraySeed
#' @return currently returns list(NULL, NULL) as we do not store dimnames in HDF5
#' @export
setMethod("dimnames", "H5S_ArraySeed", function(x) {
  list(NULL, NULL)
})
#' HDF Server content is assumed transposed relative to R matrix layout
#' @param x instance of H5S_ArraySeed
#' @return integer(2) vector of dimensions corresponding to R's layout, assuming 2-d data
#' @export
setMethod("dim", "H5S_ArraySeed", function(x) {
  # note that for HDF Server the internal dims are
  # transposed relative to R expectations
  rev(as.integer(x@H5S_dataset@shapes$dims))
})

#
# this seems incredibly convoluted 
#
setMethod("extract_array", "H5S_ArraySeed", function(x, index) {
  stopifnot(length(index)==2)
  dims_at_source = rev(dim(x))
  index = rev(index)
  full_inds_for_H5S = list(as.integer(seq_len(dims_at_source[[1]])),
        as.integer(seq_len(dims_at_source[[2]])))
  request = vector("list", 2)
# the request will operate in transpose space of R index
  if (is.null(index[[1]])) request[[2]] = full_inds_for_H5S[[1]]
    else if (length(index[[1]]) == 0) request[2] = list(NULL)
    else request[[2]] = index[[1]]
  if (is.null(index[[2]])) request[[1]] = full_inds_for_H5S[[2]]
    else if (length(index[[2]]) == 0) request[1] = list(NULL)
    else request[[1]] = index[[2]]
  reqcl = unlist(lapply(request, is.null))
  if (any(reqcl)) { # will have a zero dimension, can return without query
    if (all(reqcl)) return(matrix(numeric(0), 0, 0))
    else if (reqcl[1]) return(t(matrix(numeric(0), length(request[[2]]), 0)))
    else return(t(matrix(numeric(0), 0, length(request[[1]]))))
    }
  if (max(request[[1]]) > dims_at_source[2]) stop("dim1 request out of bounds")
  if (max(request[[2]]) > dims_at_source[1]) stop("dim2 request out of bounds")
  t(x@H5S_dataset[ request[[2]], request[[1]] ]) # native method has drop=FALSE
})
  
#' extension of DelayedArray for HDF Server content
#' @exportClass H5S_Array
setClass("H5S_Array", contains="DelayedArray")
#' extension of DelayedMatrix for HDF Server content
#' @exportClass H5S_Matrix
setClass("H5S_Matrix", contains=c("DelayedMatrix", 
     "H5S_Array"))

setMethod("matrixClass", "H5S_Array", function(x) "H5S_Matrix")

#' coercion for remote array to remote matrix

#' @name as
#' @rdname H5S_Array-class
#' @aliases coerce,H5S_Array,H5S_Matrix-method
#' @export
setAs("H5S_Array", "H5S_Matrix", function(from)
   new("H5S_Matrix", from))

setMethod("DelayedArray", "H5S_ArraySeed",
   function(seed) DelayedArray:::new_DelayedArray(seed, Class="H5S_Array"))

#' create H5S_Array instance given url (filepath) and entity (host) name
#' @param filepath a character(1) URL to port for HDF Server
#' @param host a character(1) name of 'host' in server
#' @return an instance of \code{\link[DelayedArray]{DelayedArray-class}}
#' @examples
#' # The true values from yriMulti data element 'banovichSE':
#' # > assay(banovichSE[c(1:5,329465:329469),c(1:3,63:64)])
#' #                    NA18498    NA18499    NA18501  |    NA18489     NA18909
#' # cg00000029      0.47339629  1.2943041 -0.8084735  |  0.6708168 -0.86093022
#' # cg00000165      1.23640861  0.2099817 -0.2683763  |  0.4446088  0.99868231
#' # cg00000236     -0.22258183  1.6236857 -0.8654838  |  0.1958195 -0.06090929
#' # cg00000289      0.65720581  0.5527470 -1.8458295  | -0.4618782  0.34934164
#' # cg00000363     -0.15063083  0.7498020  0.3254333  |  0.7342878  0.12940774
#' # #-------------------------------------------------------------------------
#' # ch.9.98936572R -0.07954958  0.2139431 -0.4719621  |  0.6835012  0.57758798
#' # ch.9.98937537R  0.04254705  1.0702770  1.7356387  | -0.1531732 -1.52889773
#' # ch.9.98959675F -1.59253143  0.2982456 -1.1954030  | -1.3703135  0.28974909
#' # ch.9.98989607R -1.80646652  0.4760022  1.4771808  |  0.9479602  0.49921375
#' # ch.9.991104F    0.08180195 -0.2434306  1.0281002  | -0.1653721  0.55612215
#' #
#' # compare to that delivered by H5S_Array
#' #
#' H5S_Array("http://h5s.channingremotedata.org:5000", "assays")
#' @export
H5S_Array = function(filepath, host) 
  DelayedArray(H5S_ArraySeed(filepath, host))
