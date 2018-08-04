# 
# implementation of 
# https://github.com/Bioconductor/DelayedArray/blob/master/vignettes/02-Implementing_a_backend.Rmd
# for BigQuery server back end, table with matrix layout
#
#' Represent information about a BigQuery resource with a 'triple' database schema.
#' This is targeting the isb-cgc TCGA layout.
#' BigQuery Records are regarded as triples, within major groups defined by filtervbl.
#' Triples have content subject - gene - value, to be pivoted to genes(rows) x 
#' subjects(columns) with values as entries.
#' @importFrom dplyr select_ filter_ group_by_ summarise tbl n
#' @importFrom Biobase selectSome
#' @export
setClass("BQM_Source", representation(
  bqconn = "BigQueryConnection",
  tblnm = "character",
  rowkeyfield = "character",
  allrownames = "character",
  allcolnames = "character"))
setMethod("show", "BQM_Source", function(object) {
 cat(sprintf("BQM_Source for project %s, dataset %s\n",
              object@bqconn@project, object@bqconn@dataset))
 cat("table:\n")
 cat(sprintf("\t%s\n", object@tblnm))
 cat(sprintf("rownames (of %d):\n", length(object@allrownames)))
 cat("\t", selectSome(object@allrownames), "\n")
 cat(sprintf("colnames (of %d):\n", length(object@allcolnames)))
 cat("\t", selectSome(object@allcolnames), "\n")
})

#' generate a connection to BigQuery for specific dataset
#' @note You will need to authenticate with Google.
#' @param dataset character(1) name of dataset in project
#' @param project character(1) name of project
#' @param billing character(1) billing code for project
#' @return an instance of BigQueryConnection
#' @examples
#' bqConn
#' @export
bqConn = function(dataset, project, billing) {
    requireNamespace("DBI")
    requireNamespace("bigrquery")
    con <- DBI::dbConnect(bigquery(), project = project,
        dataset = dataset, billing = billing)
    con
}

#' construct a BigQuery resource interface, for a matrix-like table
#' with one column devoted to row identification (rowkeyfield)
#' and all other columns assumed numeric
#' @param bqconn instance of BigQueryConnection from bigrquery
#' @param tblnm character(1) table name known to bqconn
#' @param rowkeyfield character(1) field in the table that will
#' @param maxdfsize numeric(1) field in the table that will constrain as.data.frame for determining rownames
#' generate rownames for matrices derived from table
#' @return instance of BQM_Source
#' @examples
#' if (interactive()) {
#'  con = bqConn(dataset="yriMulti", project=Sys.getenv("CGC_BILLING"),
#'        billing=Sys.getenv("CGC_BILLING"))
#'  banoMeth = BQM_Source(con, tblnm="banovichSE_MethylationData",
#'      rowkeyfield="cg_methyl450")
#'  banoMeth
#' }
#' @export  
BQM_Source = function(bqconn, tblnm, 
 rowkeyfield, maxdfsize=500000 ) {
 stopifnot(tblnm %in% dbListTables(bqconn))
 options(useFancyQuotes=FALSE)
 owarn = options()$warn
 options(warn=-1)
 on.exit(options(warn=owarn))
 ini = bqconn %>% tbl(tblnm) %>% select_(rowkeyfield)
 rowdf = ini %>% 
    select_(rowkeyfield) %>% group_by_(rowkeyfield) %>% summarise(n=n()) %>% as.data.frame(n=maxdfsize)
 cn = bqconn %>% tbl(tblnm) %>% as.data.frame(n=1)
 cn = setdiff(names(cn), rowkeyfield)
 new("BQM_Source", bqconn=bqconn, tblnm = tblnm,
       rowkeyfield=rowkeyfield, 
       allrownames = as.character(rowdf[,rowkeyfield]),
       allcolnames = cn)
}
#
#  
#
#' BQM_Array for BigQuery matrix content
#' @import DelayedArray
setClass("BQM_ArraySeed",
   contains="Array",
   slots = c(
     filepath="BQM_Source"))
#
BQM_ArraySeed = function(filepath) {
  requireNamespace("bigrquery")
#  tst = try(validObject(obj <- BQM_Source(filepath@bqconn)))
#  if (!is(obj, "BQM_Source")) stop("could not resolve BQM_Source request")
  stopifnot(is(filepath, "BQM_Source"))
  stopifnot(is(filepath@bqconn, "BigQueryConnection"))
  new("BQM_ArraySeed", filepath=filepath)
  }
#
#' dimnames are saved in the BQM_ArraySeed
#' @param x instance of BQM_ArraySeed
#' @return currently returns list(NULL, NULL) as we do not store dimnames in HDF5
#' @export
setMethod("dimnames", "BQM_ArraySeed", function(x) {
  list(x@filepath@allrownames, x@filepath@allcolnames)
})
#' dim derived from saved dimnames
#' @param x instance of BQM_ArraySeed
#' @return integer(2) vector of dimensions corresponding to R's layout, assuming 2-d data
#' @export
setMethod("dim", "BQM_ArraySeed", function(x) {
#  # note that for HDF Server the internal dims are
  # transposed relative to R expectations
  as.integer(c(length(x@filepath@allrownames), length(x@filepath@allcolnames)))
})
#

# this function will accept numeric or NULL for i and j
# it will generate the full data for a tumor but retrieval
# can be very slow so should be used for targeted queries
# or try to reshape in BigQuery
BQMmatgen = function(x, i, j, maxrow=Inf) {
  stopifnot(is.numeric(i) | is.null(i), is.numeric(j) | is.null(j))
  bqconn = x@filepath@bqconn
  tblnm = x@filepath@tblnm
  rowkeyfield= x@filepath@rowkeyfield
  allrows = FALSE
  allcols = FALSE
  if (!is.null(i) & length(i)>0) 
    rowsel = x@filepath@allrownames[i]
  else if (is.null(i)) {
    rowsel = x@filepath@allrownames
    allrows = TRUE  # condition the filter
    }
  else if (length(i)==0) {
    if (length(j)==0 & !is.null(j)) return(matrix(0, nrow=0, ncol=0))
      else if (is.null(j)) cn = x@filepath@allcolnames
      else cn = x@filepath@allcolnames[j]
    ans = matrix(0, nrow=0, ncol=length(cn))
    colnames(ans) = cn
    return(ans)
    }
  if (!is.null(j) & length(j)>0) {
      colsel = x@filepath@allcolnames[j]
      }
  else if (is.null(j)) {
      colsel = x@filepath@allcolnames
      allcols = TRUE
      }
  else if (length(j)==0) {
      ans = matrix(0, nrow=length(i), ncol=0)
      rownames(ans) = rowsel
      return(ans)
      }
  options(useFancyQuotes=FALSE)
  owarn = options()$warn
  options(warn=-1)
  on.exit(options(warn=owarn))
  datf = bqconn %>% tbl(tblnm) 
#
  if (!allrows) datf = datf %>%
       filter_(paste(c(rowkeyfield, "%in% rowsel"), collapse="")) # row confinement
  rn = (datf %>% select_(rowkeyfield) %>% as.data.frame(n=maxrow))[[1]]
  if (!allcols) datf = datf %>% dplyr::select(one_of(c(rowkeyfield,colsel)))
#tbl(one_of(c(rowkeyfield, colsel))) # col confinement
  datf = (datf %>% as.data.frame(n=maxrow))
  options(warn=owarn)
  rownames(datf) = as.character(datf[,rowkeyfield])
  rkind = which(colnames(datf)==rowkeyfield)
  datf = datf[,-rkind]
  mat = data.matrix(datf)
  mat[] = as.double(mat)
  mat
}


setMethod("extract_array", "BQM_ArraySeed", function(x, index) {
  stopifnot(length(index)==2)
  BQMmatgen(x, index[[1]], index[[2]], maxrow=Inf)
})
#  
#' extension of DelayedArray for BigQuery content
#' @exportClass BQM_Array
setClass("BQM_Array", contains="DelayedArray")
#' extension of DelayedMatrix for HDF Server content
#' @exportClass BQM_Matrix
setClass("BQM_Matrix", contains=c("DelayedMatrix", 
     "BQM_Array"))

setMethod("matrixClass", "BQM_Array", function(x) "BQM_Matrix")


#' coercion for remote array to remote matrix
#' @rdname BQM_Array-class
#' @aliases coerce,BQM_Array,BQM_Matrix-method
#' @import DelayedArray
#' @export
setAs("BQM_Array", "BQM_Matrix", function(from)
   new("BQM_Matrix", from))

setMethod("DelayedArray", "BQM_ArraySeed",
   function(seed) newDA(seed, Class="BQM_Array"))
#   function(seed) DelayedArray:::new_DelayedArray(seed, Class="BQM_Array"))
#
#' create BQM_Array instance given url (filepath) and entity (host) name
#' @param filepath a BQM_Source instance
#' @return an instance of \code{\link[DelayedArray]{DelayedArray-class}}
#' @examples
#'
#' # authentication issues may arise.  if you are authorized
#' # to use bigquery with GPC project isb-cgc, a token may
#' # be generated through the following
#' # options(httr_oob_default=TRUE)
#' # example(BQM_Source)
#' # a browser authentication event may occur, or if you are in
#' # a browserless session, a URL will be emitted, possibly in
#' # the context of warnings ... browse to this URL and an
#' # authentication event will occur, and a token will be provided
#' # this can be provided back to the R session to allow the
#' # query to proceed
#' #
#' if (interactive()) {
#'   con = bqConn(dataset="yriMulti", project=Sys.getenv("CGC_BILLING"),
#'        billing=Sys.getenv("CGC_BILLING"))
#'   ss = BQM_Source(con, "banovichSE_methylationData", "cg_Methyl450")
#'   #BQM_Array(ss)
#'   BQM_Array(ss)["cg00000029",c("NA18498", "NA18499", "NA18501"),drop=FALSE]
#' }
#' @export
BQM_Array = function(filepath)
  DelayedArray(BQM_ArraySeed(filepath))
