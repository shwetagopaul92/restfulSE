#' illustrate DelayedArray assay from BigQuery backend in SummarizedExperiment
#' @importFrom rlang sym UQ
#' @param bqcon a BigQueryConnection instance
#' @param colDataTableName character(1) defaulting to "clinical_PANCAN_patient_with_followup"
#' @param colDSubjectIdName character(1) defaulting to "bcr_patient_barcode"
#' @param colDFilterField character(1) defaulting to "acronym"
#' @param colDFilterValue character(1) defaulting to "BRCA"
#' @param assayDataTableName character(1) defaulting to "pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16_annot"
#' @param assayFeatureName character(1) defaulting to "ID"
#' @param assaySampleTypeCode character(1) defaulting to "TP"
#' @param subjectIDName character(1) defaulting to "ParticipantBarcode"
#' @param tumorFieldName character(1) defaulting to "Study"
#' @param tumorFieldValue character(1) defaulting to "BRCA"
#' @param assayValueFieldName character(1) defaulting to "miRNAexpr"
#' @note The parameters need different assignments for different tables.
#' Field names are not standardized across tables.
#' @return an instance of SummarizedExperiment
#' @examples
#' if (interactive() & requireNamespace("DBI") & 
#'    requireNamespace("bigrquery")) {
#'   pancan_BQ_loc = function (dataset="Annotated", 
#'                 billing=Sys.getenv("CGC_BILLING")) 
#'   {
#'       con <- DBI::dbConnect(bigrquery::bigquery(), 
#'          project = "pancancer-atlas",         
#'          dataset = dataset, billing = billing)
#'       con     
#'   }
#'   bqcon = try(pancan_BQ_loc()) # needs CGC_BILLING in environment
#'   if (!inherits(bqcon, "try-error"))
#'      pancan_SE(bqcon)
#' }
#' @export
pancan_SE = function(bqcon,
  colDataTableName = "clinical_PANCAN_patient_with_followup",
  colDSubjectIdName = "bcr_patient_barcode",
  colDFilterField = "acronym",
  colDFilterValue = "BRCA",
  assayDataTableName = "pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16_annot",
  assayFeatureName = "ID",
  assaySampleTypeCode = "TP",
  subjectIDName = "ParticipantBarcode",
  tumorFieldName = "Study",
  tumorFieldValue = "BRCA",
  assayValueFieldName = "miRNAexpr") {
 clinDF = bqcon %>% tbl(colDataTableName) %>%
   filter(UQ(rlang::sym(colDFilterField)) == colDFilterValue) %>% as.data.frame(n=Inf)
 rownames(clinDF) = clinDF[[colDSubjectIdName]]
 clinDF = DataFrame(clinDF)
 assaySrc = BQ3_Source(bqconn = bqcon,
   tblnm = assayDataTableName,
   rowkeyfield = assayFeatureName,
   assaysampletype = assaySampleTypeCode,
   colkeyfield = subjectIDName,
   filtervbl = tumorFieldName,
   filterval = tumorFieldValue,
   assayvbl = assayValueFieldName)
 assayIni = BQ3_Array(assaySrc)
 okids = intersect(rownames(clinDF), colnames(assayIni))
 se = SummarizedExperiment(assayIni[, okids], colData = clinDF[okids, ])
 names(assays(se)) = "assay"
 se
}

