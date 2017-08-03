#' simplify connection to a BigQuery dataset for ISB CGC
#' @param dataset character string with dataset name
#' @param project character string with project name
#' @param billing string with billing code
#' @export
cgcConn = function(dataset="TCGA_bioclin_v0", project="isb-cgc", billing=.cgcBilling)
  DBI::dbConnect(dbi_driver(), project = project,
        dataset = dataset, billing = billing)

#' given a BigQueryConnection to the 2016 ISB TCGA bigtables, obtain a SummarizedExperiment 'shell' rowData and colData
#' @importFrom dplyr tbl filter
#' @importFrom magrittr "%>%"
#' @param tumorCode one of the concise TCGA codes in a character string -- not checked, defaults to "LUAD", lung adenocarcinoma
#' @param assayTblName the name of the assay whose annotation will be used as rowData
#' @param rdColsToKeep columns of assay table to use in rowData component
#' @param bqConn instance of BigQueryConnection from bigrquery
#' @export
seByTumor_2016 = function(tumorCode = "LUAD", assayTblName="mRNA_UNC_HiSeq_RSEM",
    rdColsToKeep = c("original_gene_symbol", "HGNC_gene_symbol", 
      "gene_id", "Study"), bqConn ){
#
# assumes bqConn is a live instance of bigrquery BigQueryConnection
#
 stopifnot(bqConn@dataset == "tcga_201607_beta")
 require(SummarizedExperiment)
 clin = bqConn %>% tbl("Clinical_data") %>% filter(Study==tumorCode) 
 h = clin %>% select(ParticipantBarcode) %>% head() %>% as.data.frame()
 keeper = h[[1]][1]
 assay = bqConn %>% tbl(assayTblName) %>% filter(Study==tumorCode & ParticipantBarcode==keeper) %>% select(rdColsToKeep)
 SummarizedExperiment(colData=as.data.frame(clin), rowData=DataFrame(assay %>% as.data.frame()))
}

#' @exportClass BQSummarizedExperiment
setClass("BQSummarizedExperiment", contains="SummarizedExperiment",
   representation(rowQref = "ANY", colQref = "ANY",
   rowkey="character", colkey="character", assayvbl = "character"))

#' given a BigQueryConnection to the 2017 GDC-oriented ISB TCGA bigtables, obtain a SummarizedExperiment 'shell' rowData and colData
#' @param tumorCode one of the concise TCGA codes in a character string -- not checked, defaults to "LUAD", lung adenocarcinoma
#' @param assayTblName the name of the assay whose annotation will be used as rowData
#' @param rdColsToKeep columns of assay table to use in rowData component
#' @param bqConnClinical instance of BigQueryConnection from bigrquery, for access to clinical metadata -- current expectation is that the BigQuery dataset is named "TCGA_bioclin_v0" and has a table called "Clinical"
#' @param bqConnAssay instance of BigQueryConnection from bigrquery -- current expectation is that the BigQuery dataset is named "TCGA_hg19_data_v0"
#' @examples
#' require(bigrquery)
#' # be sure that .cgcBilling is set in .GlobalEnv
#' if (exists(".cgcBilling")) {
#'  clinQ = cgcConn()
#'  assayQ = cgcConn( dataset = "TCGA_hg38_data_v0" )
#'  myexpShell = seByTumor( bqConnClinical=clinQ,
#'        bqConnAssay=assayQ)
#'  print(nrow(myexpShell) == 60483)
#'  print(ncol(myexpShell) == 522)
#'  assay(myexpShell[11:15,1:4]) # some case_barcodes repeat
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
 require(SummarizedExperiment)
 newn = paste0("TCGA-", tumorCode)
 clinRef = bqConnClinical %>% tbl("Clinical")
 clin = clinRef %>% filter(project_short_name==newn)
 h = clin %>% select(case_barcode) %>% head() %>% as.data.frame()
 keeper = h[[1]][1]
 assayRef = bqConnAssay %>% tbl(assayTblName)
 assay = assayRef %>% filter(case_barcode==keeper) %>% select(rdColsToKeep)
 se = SummarizedExperiment(colData=as.data.frame(clin), rowData=DataFrame(assay %>% as.data.frame()))
 new("BQSummarizedExperiment", se, rowQref = assayRef,
     colQref = clinRef, rowkey=rowkey, colkey=colkey, assayvbl=assayvbl)
}

#' extract assay data
#' @importFrom reshape2 dcast
#' @note Very experimental approach to retrieving numerical data given
#' a SummarizedExperiment 'shell'.  We need more checking of
#' consistency between assay and clinical data before creating the
#' shell.  We use dcast to transform query result to a matrix, and
#' some 'individuals' may have multiple contributions ... we use
#' \code{fun.aggregate = max} and will see warnings until this is cleared up.
#' @exportMethod assay
setMethod("assay", c("BQSummarizedExperiment", "missing"), 
  function(x, i, ...) {
  rd = rowData(x)[,x@rowkey]
  cd = colData(x)[,x@colkey]
  df = x@rowQref %>% select(one_of(c(x@rowkey, x@colkey, x@assayvbl))) %>% filter_(paste0(x@rowkey, " %in% rd & ", x@colkey, " %in% cd")) %>% # assumes colkey shared to assay table
      as.data.frame()
  res = dcast(df, as.formula(paste0(x@rowkey, "~", x@colkey)), fun.aggregate=max, value.var=x@assayvbl)
  mat = data.matrix(res[,-1])
  rownames(mat) = as.character(res[,1])
  mat
})
