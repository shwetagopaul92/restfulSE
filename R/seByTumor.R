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

#' given a BigQueryConnection to the 2017 GDC-oriented ISB TCGA bigtables, obtain a SummarizedExperiment 'shell' rowData and colData
#' @param tumorCode one of the concise TCGA codes in a character string -- not checked, defaults to "LUAD", lung adenocarcinoma
#' @param assayTblName the name of the assay whose annotation will be used as rowData
#' @param rdColsToKeep columns of assay table to use in rowData component
#' @param bqConnClinical instance of BigQueryConnection from bigrquery, for access to clinical metadata -- current expectation is that the BigQuery dataset is named "TCGA_bioclin_v0" and has a table called "Clinical"
#' @param bqConnAssay instance of BigQueryConnection from bigrquery -- current expectation is that the BigQuery dataset is named "TCGA_hg19_data_v0"
#' @examples
#' \dontrun{
#' require(bigrquery)
#' # provide a valid billing code in the following and authenticate as needed
#' # clinQ = DBI::dbConnect(dbi_driver(), project = "isb-cgc", 
#' #       dataset = "TCGA_bioclin_v0", billing = "cgc-xx-xxxx")
#' # assayQ = DBI::dbConnect(dbi_driver(), project = "isb-cgc", 
#' #       dataset = "TCGA_hg19_data_v0", billing = "cgc-xx-xxxx")
#'  myexpShell = seByTumor_GDC( bqConnClinical=clinQ,
#'       bqConnAssay=assayQ)
#'  nrow(myexpShell) == 20531
#'  ncol(myexpShell) == 522
#' }
#' @export
seByTumor_GDC = function(tumorCode = "LUAD", 
      assayTblName="RNAseq_Gene_Expression_UNC_RSEM",
      rdColsToKeep = c("original_gene_symbol", "HGNC_gene_symbol", 
      "gene_id"), bqConnClinical, bqConnAssay ){
#
# assumes bqConn is a live instance of bigrquery BigQueryConnection
#
 if (bqConnClinical@dataset != "TCGA_bioclin_v0") warning("bqConnClinical dataset name not 'TCGA_bioclin_v0', table/column-name expectations may fail")
 if (bqConnAssay@dataset != "TCGA_hg19_data_v0") warning("bqConnClinical dataset name not 'TCGA_hg19_data_v0', table/column-name expectations may fail")
 require(SummarizedExperiment)
 newn = paste0("TCGA-", tumorCode)
 clin = bqConnClinical %>% tbl("Clinical") %>% filter(project_short_name==newn)
 h = clin %>% select(case_barcode) %>% head() %>% as.data.frame()
 keeper = h[[1]][1]
 assay = bqConnAssay %>% tbl(assayTblName) %>% filter(case_barcode==keeper) %>% select(rdColsToKeep)
 SummarizedExperiment(colData=as.data.frame(clin), rowData=DataFrame(assay %>% as.data.frame()))
}
