#' Simplify connection to a BigQuery dataset for the project "isb-cgc"
#' @importFrom bigrquery dbConnect dbi_driver
#' @param dataset character string with dataset name
#' @param project character string with project name
#' @param billing character(1) with billing code
#' @note This function operates on a BigQuery project to select a
#' dataset and return a connection.  If the google billing code is
#' assigned to environment variable \code{CGC_BILLING}, that
#' will be used to authenticate the user and collect charges.
#' Alternately the billing code can be given as a parameter.
#' @return instance of \code{\link[bigrquery:DBI]{BigQueryConnection-class}}
#' @examples
#' if (interactive()) {
#'  cgcConn
#'  require(bigrquery)
#'  # defaults concern new GDC-compliant format
#'  if (nchar(Sys.getenv("CGC_BILLING"))>0) {
#'      clin = cgcConn()
#'      dbListTables(clin)
#'  }
#' }
#' @export
cgcConn = function(dataset="TCGA_bioclin_v0", project="isb-cgc", billing=Sys.getenv("CGC_BILLING"))
  DBI::dbConnect(dbi_driver(), project = project,
        dataset = dataset, billing = billing)

#
# the following function is obsolete and is not exported
#
#' Given a BigQueryConnection to the 2016 ISB TCGA bigtables, obtain a SummarizedExperiment 'shell' rowData and colData
#' @importFrom dplyr tbl filter filter_ one_of
#' @importFrom magrittr "%>%"
#' @import SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom utils head
#' @param tumorCode one of the concise TCGA codes in a character string -- not checked, defaults to "LUAD", lung adenocarcinoma
#' @param assayTblName the name of the assay whose annotation will be used as rowData
#' @param rdColsToKeep columns of assay table to use in rowData component
#' @param bqConn instance of BigQueryConnection from bigrquery
#' @return SummarizedExperiment instance, with BigQuery reference as assay
seByTumor_2016 = function(tumorCode = "LUAD", assayTblName="mRNA_UNC_HiSeq_RSEM",
    rdColsToKeep = c("original_gene_symbol", "HGNC_gene_symbol", 
      "gene_id", "Study"), bqConn ){
#
# assumes bqConn is a live instance of bigrquery BigQueryConnection
#
 stopifnot(bqConn@dataset == "tcga_201607_beta")
 requireNamespace("SummarizedExperiment")
 Study = NULL
 ParticipantBarcode = NULL
 clin = bqConn %>% tbl("Clinical_data") %>% filter(Study==tumorCode) 
 h = clin %>% dplyr::select(ParticipantBarcode) %>% head() %>% as.data.frame()
 keeper = h[[1]][1]
 assay = bqConn %>% tbl(assayTblName) %>% filter(Study==tumorCode & ParticipantBarcode==keeper) %>% dplyr::select(rdColsToKeep)
 SummarizedExperiment(colData=as.data.frame(clin), rowData=DataFrame(assay %>% as.data.frame()))
}

#' Define a class to use BigQuery data through SummarizedExperiment interface
#' @rdname BQSummarizedExperiment
#' @slot rowQref a BigQueryConnection wrapped in tbl_dbi that
#' holds rowData for the SummarizedExperiment instance
#' @slot colQref a BigQueryConnection wrapped in tbl_dbi that
#' holds colData for the SummarizedExperiment instance
#' @slot rowkey character(1) name of a field in the table
#' referenced by \code{rowQref} to be used as key for features
#' @slot colkey character(1) name of a field in the table
#' referenced by \code{colQref} to use as key for samples
#' @slot assayvbl character(1) name to be used to select table
#' providing assay content
#' @note This is an experimental structure to probe the concept that
#' one can use a SummarizedExperiment object to interact with
#' BigQuery data, particularly TCGA data.  The slots \code{rowQref}
#' and \code{colQref} are expected to be BigQuery connections
#' which supply information on features and samples respectively, in
#' a way that is consistent with the assay representation.
#' See \code{\link{seByTumor}} for illustration.
#' @exportClass BQSummarizedExperiment
setClass("BQSummarizedExperiment", contains="SummarizedExperiment",
   representation(rowQref = "ANY", colQref = "ANY",
   rowkey="character", colkey="character", assayvbl = "character"))

#' Given a BigQueryConnection to the 2017 GDC-oriented ISB TCGA bigtables, obtain a SummarizedExperiment 'shell' rowData and colData
#' @param tumorCode one of the concise TCGA codes in a character string -- not checked, defaults to "LUAD", lung adenocarcinoma
#' @param assayTblName the name of the assay whose annotation will be used as rowData
#' @param rdColsToKeep columns of assay table to use in rowData component
#' @param bqConnClinical instance of BigQueryConnection from bigrquery, for access to clinical metadata -- current expectation is that the BigQuery dataset is named "TCGA_bioclin_v0" and has a table called "Clinical"
#' @param bqConnAssay instance of BigQueryConnection from bigrquery -- current expectation is that the BigQuery dataset is named "TCGA_hg19_data_v0"
#' @param rowkey name of a field to be used as key for rows
#' @param colkey name of a field to use as key for samples
#' @param assayvbl name of field to use for numerical values
#' @return SummarizedExperiment
#' @note This function demonstrates the use of external resources
#' for rowData, colData and assay components of a SummarizedExperiment
#' instance.  The intention is that the full complement of activities
#' supported by \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
#' are likewise supported through this class, with 
#' assay data and sample and feature
#' metadata all external and in BigQuery projects.  The \code{seByTumor}
#' function is provided to generate an example of this approach
#' with minimal user configuration.
#' @examples
#' if (interactive()) {
#'  require(bigrquery)
#'  # be sure that .cgcBilling is set 
#'  code = Sys.getenv("CGC_BILLING")
#'  if (!(nchar(code)==0)) {
#'   clinQ = cgcConn(billing=code)
#'   assayQ = cgcConn( dataset = "TCGA_hg38_data_v0", billing=code )
#'   myexpShell = seByTumor( bqConnClinical=clinQ,
#'         bqConnAssay=assayQ)
#'   print(myexpShell)
#'   print(nrow(myexpShell) == 60483)
#'   print(ncol(myexpShell) == 522)
#'   assay(myexpShell[11:15,1:4]) # some case_barcodes repeat
#'   }
#'  }
#' @export
seByTumor = function(tumorCode = "LUAD", 
      assayTblName="RNAseq_Gene_Expression",
      rdColsToKeep = c("gene_name", "Ensembl_gene_id", "gene_type"),
      bqConnClinical, bqConnAssay, rowkey="Ensembl_gene_id",
      colkey = "case_barcode", assayvbl="HTSeq__Counts" ){
#
# assumes bqConn is a live instance of bigrquery BigQueryConnection
#
 if (bqConnClinical@dataset != "TCGA_bioclin_v0") warning("bqConnClinical dataset name not 'TCGA_bioclin_v0', table/column-name expectations may fail")
 if (bqConnAssay@dataset != "TCGA_hg38_data_v0") warning("bqConnAssay dataset name not 'TCGA_hg38data_v0', table/column-name expectations may fail")
 case_barcode = NULL
 project_short_name = NULL

# suggestion from mike lawrence
#clin <- clinRef %>% filter(project_short_name==newn)
#rowData <- inner_join(clin, assayRef, by=colkey) %>%
#    dplyr::select(rdColsToKeep)

 requireNamespace("SummarizedExperiment")
 newn = paste0("TCGA-", tumorCode)
 clinRef = bqConnClinical %>% tbl("Clinical")
 clin = clinRef %>% filter(project_short_name==newn)
 h = clin %>% dplyr::select(case_barcode) %>% head() %>% as.data.frame()
 keeper = h[[1]][1]
 assayRef = bqConnAssay %>% tbl(assayTblName)
 assay = assayRef %>% filter(case_barcode==keeper) %>% dplyr::select(rdColsToKeep)
 clin = as.data.frame(clin)
 rownames(clin) = clin[,colkey]
 se = SummarizedExperiment(colData=clin, rowData=DataFrame(assay %>% as.data.frame()))
 se@NAMES = rowData(se)[,rowkey]
 new("BQSummarizedExperiment", se, rowQref = assayRef,
     colQref = clinRef, rowkey=rowkey, colkey=colkey, assayvbl=assayvbl)
}

#' extract assay data
#' @importFrom reshape2 dcast
#' @importFrom stats as.formula
#' @param x BQSummarizedExperiment instance
#' @param i index for retrieval, ignored at present
#' @param \dots not used
#' @note Very experimental approach to retrieving numerical data given
#' a SummarizedExperiment 'shell'.  We need more checking of
#' consistency between assay and clinical data before creating the
#' shell.  We use dcast to transform query result to a matrix, and
#' some 'individuals' may have multiple contributions ... we use
#' \code{fun.aggregate = max} and will see warnings until this is cleared up.
#' @return matrix
#' @exportMethod assay
setMethod("assay", c("BQSummarizedExperiment", "missing"), 
  function(x, i, ...) {
  rd = rowData(x)[,x@rowkey]
  cd = colData(x)[,x@colkey]

  df = x@rowQref %>% dplyr::select(one_of(c(x@rowkey, x@colkey, x@assayvbl))) %>% filter_(paste0(x@rowkey, " %in% rd & ", x@colkey, " %in% cd")) %>% # assumes colkey shared to assay table
      as.data.frame()
  res = dcast(df, as.formula(paste0(x@rowkey, "~", x@colkey)), fun.aggregate=max, value.var=x@assayvbl)
  mat = data.matrix(res[,-1])
  rownames(mat) = as.character(res[,1])
  mat
})

#' Placeholder for assay name extractor for a BQSummarizedExperiment instance.
#' @note This function supplies a placeholder for this early
#' version of a SummarizedExperiment instance to BigQuery.
#' At present there is only one assay supported; future work will
#' help to reduce special coding for BigQuery back end.
#' @aliases assayNames,BQSummarizedExperiment-method
#' @aliases assayNames
#' @param x instance of BQSummarizedExperiment
#' @param \dots not used
#' @return string indicating that assay is served by BigQuery, nameless
#' @exportMethod assayNames
setMethod("assayNames", "BQSummarizedExperiment", function(x, ...) {
 "(served by BigQuery)"
})


