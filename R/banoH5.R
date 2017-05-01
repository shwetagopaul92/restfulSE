#' interact with banovichSE content via hdf5 server
#' @export
banoH5 = function(url= "http://54.163.220.201:5000",
          host="assays.hdfgroup.org", ...) {
#  "http://170.223.248.164:7248"
con = getH5ShContent(url, host)
ds = datasets(con)
ds
}

