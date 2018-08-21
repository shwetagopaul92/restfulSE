
doit_TCGA = function() {
library(restfulSE)
library(bigrquery)
library(DBI)
library(testthat)
library(reshape2)
 NSAMP = 3
 NFEAT = 2
 if (nchar(Sys.getenv("CGC_BILLING"))==0) skip("CGC_BILLING not set")
 con = try(cgcConn("TCGA_hg38_data_v0"))
 if (inherits(con, "try-error")) skip("can't establish CGC_BILLING")
 ss = BQ3_Source(con, filterval="TCGA-LUAD")
 samps = ss@allcolnames[sample(seq_len(length(ss@allcolnames)),NSAMP)]
 feats = ss@allrownames[sample(seq_len(length(ss@allrownames)),NFEAT)]
 myq = "SELECT ensembl_gene_id, case_barcode, HTSeq__Counts FROM `isb-cgc.TCGA_hg38_data_v0.RNAseq_Gene_Expression` where \n(case_barcode = '%s' or case_barcode = '%s' or case_barcode = '%s') and \n  project_short_name = 'TCGA-LUAD' and (ensembl_gene_id = '%s' or ensembl_gene_id = '%s')" 
 myq = sprintf(myq, samps[1], samps[2], samps[3], feats[1], feats[2])
 goodTable = DBI::dbGetQuery(con, myq)
 expect_true(nrow(goodTable)==6) # add more tests for content of
# BQ3_Array execution
 arr = BQ3_Array(ss)[c(feats[1], feats[2]),
    c(samps[1], samps[2], samps[3])]
 expect_true(inherits(arr, "DelayedArray"))
 mat = as.matrix(arr)
 goodt = dcast(ensembl_gene_id~case_barcode, data=goodTable)
 rownames(goodt) = goodt[,1]
 all.equal(data.matrix(goodt[rownames(mat), colnames(mat)]), mat)
}

#for(i in 1:20) {
#  x = doit()
#  cat(i)
#  stopifnot(x==TRUE)
#}

context("bigquery accuracy for project isb-cgc")
test_that("random query resolves directly and via BQ3_Array", {
 expect_true(doit_TCGA())
#
# note that an example with project pancan-atlas is present in BiocOncoTK
#
})
