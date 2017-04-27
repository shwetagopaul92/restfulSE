#' interact with banovichSE content via hdf5 server
#' @export
banoH5 = function(url="http://170.223.248.164:7248",
          host="assays.hdfgroup.org", ...) {
con = getH5ShContent(url, host)
ds = datasets(con)
ds
}

