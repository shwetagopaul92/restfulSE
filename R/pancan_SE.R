#' illustrate DelayedArray assay from BigQuery backend in SummarizedExperiment
#' @importFrom rlang sym UQ
#' @param bqcon a BigQueryConnection instance
#' @param colDataTableName character(1) defaulting to "clinical_PANCAN_patient_with_followup"
#' @param clinVars character() vector of variables to be retained from the table named by 'colDataTableName', defaults to vector returned by clinVars()
#' @param colDSubjectIdName character(1) defaulting to "bcr_patient_barcode"
#' @param colDFilterField character(1) defaulting to "acronym"
#' @param colDFilterValue character(1) defaulting to "BRCA"; a vector may be
#' used, in which case multiple tumor types will be 
#' represented -- must agree with tumorFieldValue, as clinical
#' and assay data are collected separately
#' @param assayDataTableName character(1) defaulting to "pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16_annot"
#' @param assayFeatureName character(1) defaulting to "ID"
#' @param assaySampleTypeCode character(1) defaulting to "TP"
#' @param subjectIDName character(1) defaulting to "ParticipantBarcode"
#' @param tumorFieldName character(1) defaulting to "Study"
#' @param tumorFieldValue character() defaulting to "BRCA"; a vector may be used in which case multiple tumor types will be represented -- must agree with colDFilterValue
#' @param assayValueFieldName character(1) defaulting to "miRNAexpr"
#' @note The parameters need different assignments for different tables.
#' Field names are not standardized across tables as of August 2018.  AUTHENTICATION CONCERNS:
#' You must have a valid BigQuery project identifier in the environment variable
#' CGC_BILLING, or pass such as 'billing' when using DBI::dbConnect(bigquery::bigrquery(), ...).
#' To get such a project identifier, you need to have a Google identity and you must
#' have created a BigQuery project with that identity.  Notes at 
#' \url{https://isb-cancer-genomics-cloud.readthedocs.io/en/latest/sections/progapi/bigqueryGUI/WalkthroughOfGoogleBigQuery.html} provide details.
#' Another approach that can be used involves the Google Cloud SDK.  Assuming you have
#' a Google identity and have installed a recent SDK instance, you can use (in the shell)
#' \code{gcloud auth login} to pick the identity that has a project with id
#' \code{PROJECT_ID}.  Use that project id as the billing code for bigrquery dbConnect,
#' or set CGC_BILLING in the environment to evaluate to that project id.
#' @return an instance of SummarizedExperiment
#' @examples
#' if (interactive() & requireNamespace("DBI") & 
#'    requireNamespace("bigrquery")) {
#'      billco = Sys.getenv("CGC_BILLING")
#'      if (nchar(billco)==0) 
#'           stop("need CGC_BILLING set to your BigQuery project ID, see note in ?pancan_SE")
#'      bqcon = DBI::dbConnect(bigrquery::bigquery(), project = "pancancer-atlas", 
#'            dataset = "Annotated", billing = billco)
#'      brca_mirSE = pancan_SE(bqcon)
#'      brca_mirSE
#'      }
#' @export
pancan_SE = function(bqcon,
  colDataTableName = "clinical_PANCAN_patient_with_followup",
  clinVars = basic_clinvars(),
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
 stopifnot(all.equal(colDFilterValue, tumorFieldValue))
 msg = "Cannot probe requested clinical table.  This is probably an authentication
issue.  Please read the note in ?pancan_SE to establish a valid .httr-oauth
before retrying."
 clinDF = tryCatch( bqcon %>% tbl(colDataTableName) %>% select(UQ(clinVars)) %>%
   filter(UQ(rlang::sym(colDFilterField)) %in% colDFilterValue) %>% as.data.frame(n=Inf),
     error = function(e) e )
 if (inherits(clinDF, "simpleError")) stop(msg) # probably 404
 if (inherits(clinDF, "bigrquery_notFound")) {
   print(clinDF)
   stop("Probably the value of environment variable CGC_BILLING 
in the BigQueryConnection does not correspond to a valid BigQuery 
project to which pancancer-atlas tables are linked.")
   }
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

#' provide a collection of basic clinical variables to limit size of clinical data return
#' @return character(20) vector of variable names
#' @note Use pancan_app() to determine other variables available
#' @export
basic_clinvars = function() {
c("bcr_patient_uuid", "bcr_patient_barcode", "acronym", "gender", 
"vital_status", "days_to_birth", "days_to_death", "days_to_last_followup", 
"days_to_initial_pathologic_diagnosis", "age_at_initial_pathologic_diagnosis", 
"icd_10", "tissue_retrospective_collection_indicator", "icd_o_3_histology", 
"tissue_prospective_collection_indicator", "history_of_neoadjuvant_treatment", 
"icd_o_3_site", "tumor_tissue_site", "new_tumor_event_after_initial_treatment", 
"radiation_therapy", "race")
}
