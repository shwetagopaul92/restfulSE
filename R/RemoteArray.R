# longrunning jan 2018
# domain http://54.174.163.77:5000
# host tenx_full

#' RemoteArray for HDF Server content
#' @import DelayedArray
#' @exportClass RemoteArray
#' @exportClass RemoteMatrix
setClass("RemoteArraySeed",
   contains="Array",
   slots = c(
     filepath="character",
     domain="character",
     host="character",
     H5S_dataset="H5S_dataset"
     ))

RemoteArraySeed = function(filepath, host) {
  requireNamespace("rhdf5client")
  dom = try(rhdf5client::H5S_source(filepath))
  if (!is(dom, "H5S_source")) stop("could not resolve H5S_source request on filepath")
  ds = try(dom[[host]])
  if (!is(ds, "H5S_dataset")) stop("could not resolve H5S_dataset request on filepath[[host]]")
  new("RemoteArraySeed", filepath=filepath,
         domain=filepath, host=host, H5S_dataset=ds)
  }

setMethod("dimnames", "RemoteArraySeed", function(x) {
  list(NULL, NULL)
})
setMethod("dim", "RemoteArraySeed", function(x) {
  # note that for HDF Server the internal dims are
  # transposed relative to R expectations
  rev(as.integer(x@H5S_dataset@shapes$dims))
})

#
# this seems incredibly convoluted -- VJ Carey, the author
#
setMethod("extract_array", "RemoteArraySeed", function(x, index) {
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
  
setClass("RemoteArray", contains="DelayedArray")
setClass("RemoteMatrix", contains=c("DelayedMatrix", 
     "RemoteArray"))

setMethod("matrixClass", "RemoteArray", function(x) "RemoteMatrix")
#' coercion for remote array to remote matrix
#' @name setAs-RemoteArray-RemoteMatrix
#' @rdname RemoteArray-class
#' @exportMethod coerce
setAs("RemoteArray", "RemoteMatrix", function(from)
   new("RemoteMatrix", from))

setMethod("DelayedArray", "RemoteArraySeed",
   function(seed) DelayedArray:::new_DelayedArray(seed, Class="RemoteArray"))

#' create RemoteArray instance given url (filepath) and entity (host) name
#' @param filepath a character(1) URL to port for HDF Server
#' @param host a character(1) name of 'host' in server
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
#' # compare to that delivered by RemoteArray
#' #
#' RemoteArray("http://h5s.channingremotedata.org:5000", "assays")
#' @export
RemoteArray = function(filepath, host) 
  DelayedArray(RemoteArraySeed(filepath, host))


  
  
