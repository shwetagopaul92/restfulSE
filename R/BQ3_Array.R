#' @rawNamespace importClassesFrom("DelayedArray", "DelayedArray")
#' @importFrom DelayedArray matrixClass
#' @importFrom S4Vectors new2
newDA = function (seed = new("array"), Class = "DelayedArray") 
{
    seed_ndim <- length(dim(seed))
    if (seed_ndim == 2L) 
        Class <- matrixClass(new(Class))
    new2(Class, seed = seed)
}

# 
# implementation of 
# https://github.com/Bioconductor/DelayedArray/blob/master/vignettes/02-Implementing_a_backend.Rmd
# for BigQuery server back end, table with 'triples'
#
#' Represent information about a BigQuery resource with a 'triple' database schema.
#' This is targeting the isb-cgc TCGA layout.
#' BigQuery Records are regarded as triples, within major groups defined by filtervbl.
#' Triples have content subject - gene - value, to be pivoted to genes(rows) x 
#' subjects(columns) with values as entries.
#' @importFrom dplyr select_ filter_ group_by_ summarise tbl n
#' @importFrom Biobase selectSome
#' @import DelayedArray
#' @export
setClass("BQ3_Source", representation(
  bqconn = "BigQueryConnection",
  tblnm = "character",
  rowkeyfield = "character",
  colkeyfield = "character",
  filterfield = "character",
  filterval = "character",
  filtervbl = "character",
  assayvbl = "character",
  assaysampletype = "character",
  allrownames = "character",
  allcolnames = "character"))
setMethod("show", "BQ3_Source", function(object) {
 cat(sprintf("BQ3_Source for project %s, dataset %s, assayvbl %s\n",
              object@bqconn@project, object@bqconn@dataset, object@assayvbl))
 cat("table [filter]:\n")
 cat(sprintf("\t%s, [%s == %s]\n", object@tblnm, object@filtervbl, sQuote(object@filterval)))
 cat(sprintf("rownames (of %d):\n", length(object@allrownames)))
 cat("\t", selectSome(object@allrownames), "\n")
 cat(sprintf("colnames (of %d):\n", length(object@allcolnames)))
 cat("\t", selectSome(object@allcolnames), "\n")
})

#' construct a BigQuery resource interface
#' @param bqconn instance of BigQueryConnection from bigrquery
#' @param tblnm character(1) table name known to bqconn
#' @param rowkeyfield character(1) field in the table that will
#' generate rownames for matrices derived from table
#' @param colkeyfield character(1) field in the table that will
#' generate colnames for matrices derived from table
#' @param filtervbl character(1) field in the table that will be used to filter out a group of records,
#' for example, all records pertaining to a given tumor in TCGA
#' @param filterval character(1) value in the range of filtervbl to identify records to retain --
#' @param assayvbl character(1) field with assay quantifications
#' @param assaysampletype character(1) value for filtering pancancer-atlas assays, which include normals and other sample types, defaulting to "TP"; ignored if project element of \code{bqconn} is not `pancancer-atlas`
#' @return instance of BQ3_Source
#' @examples
#' if (interactive()) {
#'  con = cgcConn("TCGA_hg38_data_v0")
#'  lungConn = BQ3_Source(con, filterval="TCGA-LUAD")
#'  lungConn
#' }
#' @export  
BQ3_Source = function(bqconn, tblnm = "RNAseq_Gene_Expression",
 rowkeyfield = "Ensembl_gene_id", colkeyfield = "case_barcode",
 filtervbl = "project_short_name", filterval = "TCGA-GBM",
   assayvbl = "HTSeq__Counts", assaysampletype="TP") {
 stopifnot(tblnm %in% dbListTables(bqconn))
 options(useFancyQuotes=FALSE)
 if (slot(bqconn, "project") != "pancancer-atlas") {
 ini = bqconn %>% tbl(tblnm) %>% select_(rowkeyfield, filtervbl, 
        colkeyfield) %>%
    filter_(paste(c(filtervbl, "==", sQuote(filterval)), collapse="")) 
    } else ini = 
    bqconn %>% tbl(tblnm) %>% select_(rowkeyfield, filtervbl, 
        colkeyfield, "SampleTypeLetterCode") %>%
    filter_(paste(c(filtervbl, "==", sQuote(filterval)), collapse="")) %>%
       filter(SampleTypeLetterCode == assaysampletype)
 rowdf = ini %>% 
    select_(rowkeyfield) %>% group_by_(rowkeyfield) %>% summarise(n=n()) %>% as.data.frame(n=100000)
 coldf = ini %>%
    select_(colkeyfield) %>% group_by_(colkeyfield) %>% summarise(n=n()) %>% as.data.frame(n=100000)
 colns = coldf[,2]
 ntab = table(colns)
 modal = ntab[which.max(ntab)]
 outl = which(ntab != modal)
 if (length(outl)>0) {
   message(paste(colkeyfield, "has", sum(ntab[outl]), "contributors with excess contributions that are omitted"))
   coldf = coldf[ which(coldf[,2] == as.numeric(names(modal))), ]
   }
 new("BQ3_Source", bqconn=bqconn, tblnm = tblnm,
       rowkeyfield=rowkeyfield, colkeyfield=colkeyfield,
       filtervbl = filtervbl, filterval = filterval,
       assayvbl = assayvbl,
       assaysampletype = assaysampletype,
       allrownames = as.character(rowdf[,rowkeyfield]),
       allcolnames = as.character(coldf[,colkeyfield]))
}
#
#  
#
#' BQ3_Array for BigQuery matrix content
#' @import DelayedArray
setClass("BQ3_ArraySeed",
   contains="Array",
   slots = c(
     filepath="BQ3_Source"))
#'@import DelayedArray
BQ3_ArraySeed = function(filepath) {
  requireNamespace("bigrquery")
#  tst = try(validObject(obj <- BQ3_Source(filepath@bqconn)))
#  if (!is(obj, "BQ3_Source")) stop("could not resolve BQ3_Source request")
  stopifnot(is(filepath, "BQ3_Source"))
  stopifnot(is(filepath@bqconn, "BigQueryConnection"))
  new("BQ3_ArraySeed", filepath=filepath)
  }
#
#' dimnames are saved in the BQ3_ArraySeed
#' @param x instance of BQ3_ArraySeed
#' @return currently returns list(NULL, NULL) as we do not store dimnames in HDF5
#' @export
setMethod("dimnames", "BQ3_ArraySeed", function(x) {
  list(x@filepath@allrownames, x@filepath@allcolnames)
})
#' dim derived from saved dimnames
#' @param x instance of BQ3_ArraySeed
#' @return integer(2) vector of dimensions corresponding to R's layout, assuming 2-d data
#' @export
setMethod("dim", "BQ3_ArraySeed", function(x) {
#  # note that for HDF Server the internal dims are
  # transposed relative to R expectations
  as.integer(c(length(x@filepath@allrownames), length(x@filepath@allcolnames)))
})
#

# this function will accept numeric or NULL for i and j
# it will generate the full data for a tumor but retrieval
# can be very slow so should be used for targeted queries
# or try to reshape in BigQuery
BQ3matgen = function(x, i, j, maxrow=Inf) {
  stopifnot(is.numeric(i) | is.null(i), is.numeric(j) | is.null(j))
  bqconn = x@filepath@bqconn
  tblnm = x@filepath@tblnm
  rowkeyfield = x@filepath@rowkeyfield
  colkeyfield = x@filepath@colkeyfield
  filtervbl = x@filepath@filtervbl
  filterval = x@filepath@filterval
  assayvbl = x@filepath@assayvbl
  assaysampletype = x@filepath@assaysampletype
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
  df = bqconn %>% tbl(tblnm) %>%   # QUESTION: Can this reference be carried in the seed?
       select_(rowkeyfield, colkeyfield, filtervbl, assayvbl) %>%  # confine columns
       filter_(paste(c(filtervbl, "==", sQuote(filterval)), collapse="")) # major row confinement
  if (!allcols) df = df %>%
       filter_(paste(c(colkeyfield, "%in% colsel"), collapse="")) # col confinement
  if (!allrows) df = df %>%
       filter_(paste(c(rowkeyfield, "%in% rowsel"), collapse="")) # minor row confinement
  df = (df %>% as.data.frame(n=maxrow))
  df = df[ which(df[[colkeyfield]] %in% x@filepath@allcolnames), ]
  df = dcast(df, as.formula(paste(rowkeyfield, "~", colkeyfield, collapse="")), value.var=assayvbl, fun.aggregate=mean)
  rownames(df) = df[,1]
  df = df[,-1]
  mat = data.matrix(df)
  mat[] = as.double(mat)
  mat
}
#' @import DelayedArray
setMethod("extract_array", "BQ3_ArraySeed", function(x, index) {
  stopifnot(length(index)==2)
  BQ3matgen(x, index[[1]], index[[2]], maxrow=Inf)
})
#  
#' extension of DelayedArray for BigQuery content
#' @exportClass BQ3_Array
setClass("BQ3_Array", contains="DelayedArray")
#' extension of DelayedMatrix for HDF Server content
#' @exportClass BQ3_Matrix
setClass("BQ3_Matrix", contains=c("DelayedMatrix", 
     "BQ3_Array"))

setMethod("matrixClass", "BQ3_Array", function(x) "BQ3_Matrix")


#' coercion for remote array to remote matrix
#' @rdname BQ3_Array-class
#' @aliases coerce,BQ3_Array,BQ3_Matrix-method
#' @import DelayedArray
#' @export
setAs("BQ3_Array", "BQ3_Matrix", function(from)
   new("BQ3_Matrix", from))

setMethod("DelayedArray", "BQ3_ArraySeed",
#   function(seed) DelayedArray:::new_DelayedArray(seed, Class="BQ3_Array"))
   function(seed) newDA(seed, Class="BQ3_Array"))
#
#' create BQ3_Array instance given url (filepath) and entity (host) name
#' @param filepath a BQ3_Source instance
#' @return an instance of \code{\link[DelayedArray]{DelayedArray-class}}
#' @examples
#'
#' # authentication issues may arise.  if you are authorized
#' # to use bigquery with GPC project isb-cgc, a token may
#' # be generated through the following
#' # options(httr_oob_default=TRUE)
#' # example(BQ3_Source)
#' # a browser authentication event may occur, or if you are in
#' # a browserless session, a URL will be emitted, possibly in
#' # the context of warnings ... browse to this URL and an
#' # authentication event will occur, and a token will be provided
#' # this can be provided back to the R session to allow the
#' # query to proceed
#' #
#' if (interactive()) {
#'   con = cgcConn("TCGA_hg38_data_v0")
#'   ss = BQ3_Source(con, filterval="TCGA-LUAD")
#'   BQ3_Array(ss)
#' }
#' @export
BQ3_Array = function(filepath)
  DelayedArray(BQ3_ArraySeed(filepath))
