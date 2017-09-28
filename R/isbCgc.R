# some utilities to simplify exploration of ISB Cancer Genomic Cloud BigQuery
# August 2017

#' vector of dataset names in isb-cgc project
#' @return character vector
#' @examples
#' isbCgcDatasets() # might be out of date ... can use list_datasets in bigrquery
#' @export
isbCgcDatasets = function() {
 c("ccle_201602_alpha",
 "GDC_metadata",
 "genome_reference",
 "hg19_data_previews",
 "hg38_data_previews",
 "metadata",
 "platform_reference",
 "QotM",
 "TARGET_bioclin_v0",
 "TARGET_hg38_data_v0",
 "tcga_201607_beta",
 "TCGA_bioclin_v0",
 "tcga_cohorts",
 "TCGA_hg19_data_v0",
 "TCGA_hg38_data_v0",
 "tcga_seq_metadata",
 "Toil_recompute")
}

#' list the tables in a selected dataset
#' @import DBI
#' @param dataset character string identifying a table in ISB CGC
#' @param billing Google BigQuery billing code, which can be set in an environment variable \code{CGC_BILLING}
#' @return character vector
#' @examples 
#' # be sure that .cgcBilling is set
#' code = Sys.getenv("CGC_BILLING")
#' if (!(nchar(code)==0)) {
#'  isbCgcTables()
#'  }
#' @export
isbCgcTables = function(dataset="TCGA_hg19_data_v0", billing=Sys.getenv("CGC_BILLING")) {
  stopifnot(dataset %in% isbCgcDatasets())
  con <- dbConnect(dbi_driver(), project = "isb-cgc", 
        dataset = dataset, billing = billing)
  on.exit(dbDisconnect(con))
  dbListTables(con)
}  

setClass("TableSet", representation(dataset="character", 
    tablenames="character", tables="list"))
setMethod("show", "TableSet", function(object) {
 cat("TableSet for dataset ", object@dataset, "\n")
 cat(" with", length(object@tables), "tables.\n")

}) 

TCGA_tablerefs = function(build="hg19", billing=Sys.getenv("CGC_BILLING")) {
 getConn = function(dataset) dbConnect(dbi_driver(), project = "isb-cgc", 
        dataset = dataset, billing = billing)
 if (build == "hg19") {
   ds = "TCGA_hg19_data_v0"
   }
 else if (build == "hg38") {
   ds = "TCGA_hg38_data_v0"
   }
 tblnames = isbCgcTables(dataset=ds, billing=billing)
 new("TableSet", dataset=ds, tablenames=tblnames, tables=
      lapply(tblnames, function(x) getConn(ds) %>% tbl(x)))
}

  
