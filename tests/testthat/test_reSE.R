library(restfulSE)

context("connection")

test_that("H5S_source completes", {
 bigec2 = H5S_source("http://54.174.163.77:5000")
 expect_true(is(bigec2, "H5S_source"))
})

context("content wrapper structure") 

test_that("H5S_source processes", {
 bigec2 = H5S_source("http://54.174.163.77:5000")
 expect_true(all(dim(groups(bigec2))==c(10,2))) 
 expect_true(is(links(bigec2,1), "H5S_linkset"))
 expect_true(is(dataset(bigec2, "tenx_100k"), "H5S_dataset"))
 expect_true(is(bigec2[["tenx_100k"]], "H5S_dataset"))
 expect_true(is(dsmeta(bigec2), "DataFrame"))
 expect_true(all(dim(dsmeta(bigec2))==c(10,3))) 
})

context("indexing infrastructure")

test_that("sproc/isplit work", {
 expect_true(length(isplit(c(1,2,3,4,5,10,15,20,30:40)))==4)
 ii = isplit(c(1,2,3,4,5,10,15,20,30:40))
 ss = structure(c("0:5:1", "9:20:5", "29:30:1", "30:40:1"), .Names = c("1", 
"2", "3", "4"))
 expect_true(identical(ss, unlist(sproc(ii))))
})

context("targets generation")

test_that("targets method works", {
 bigec2 = H5S_source("http://54.174.163.77:5000")
 tt = targets(links(bigec2, 1))
 expect_true(length(tt)==9)  # july 5 2017
 expect_true(length(grep("host", tt))==7) # added tenx_full dataset
})

context("RESTful SE")

test_that("RESTfulSummarizedExperiment infrastructure works against server", {
 bigec2 = H5S_source("http://54.174.163.77:5000")
 data(st100k)
 n100k_served = bigec2[["tenx_100k_sorted"]]  
 rr = RESTfulSummarizedExperiment( st100k, n100k_served )
 expect_true(is(rr, "RangedSummarizedExperiment"))
 data(banoSEMeta) 
 bigec2 = H5S_source("http://54.174.163.77:5000")
 banoh5 = bigec2[["assays"]]
 rr = RESTfulSummarizedExperiment( banoSEMeta, banoh5 )
 expect_true(validObject(rr))
 expect_true(all(dim(rr)==rev(internalDim(banoh5))))
 rr2 = rr[5:12000,3:10]
 arr2 = assay(rr2)
 arr2s = structure(c(310.234826991, -574.26777552, 
            -443.279313398, -8.31305800099999, 
            -420.419219739, -41.648322087, 
            160.146168073, -39.986112183), 
         .Names = c("NA18501", "NA18502", "NA18516", "NA18517", 
                     "NA18519", "NA18520", "NA18855", "NA18856"))
 csums = apply(arr2,2,sum)
 expect_true(max(abs(arr2s-csums)) < 1e-6)
# sumstr = structure(c(3667, 1844, 4329, 2928, 7787, 10084, 2151, 3318, 
#2851, 4158), .Names = c("AAACCTGAGATAGGAG-1", "AAACCTGAGCGGCTTC-1", 
#"AAACCTGAGGAATCGC-1", "AAACCTGAGGACACCA-1", "AAACCTGAGGCCCGTT-1", 
#"AAACCTGAGTCCGGTC-1", "AAACCTGCAACACGCC-1", "AAACCTGCACAGCGTC-1", 
#"AAACCTGCAGCCACCA-1", "AAACCTGCAGGATTGG-1"))
# expect_true(identical(apply(arr2,2,sum), sumstr))
})
 
test_that("complex indexing succeeds", {
 library(yriMulti)
 data(banovichSE)
 bansel = assay(banovichSE[c(1,3,5,200000,300000),c(1,5,7)])
 library(restfulSE)
 data(banoSEMeta) 
 bigec2 = H5S_source("http://54.174.163.77:5000")
 banh = bigec2[["assays"]]
 rr = RESTfulSummarizedExperiment( banoSEMeta, banh )
 rrsel = rr[c(1,3,5,200000,300000),c(1,5,7)]
 rrass = assay(rrsel)
 expect_true( max(abs(rrass-bansel))<1e-6 )
 rrsel0 = rr[c(1,3,5,200000,300000),c(1)]
 rrass0 = assay(rrsel0)
 expect_true(ncol(rrass0)==1)
})

test_that("dim compatibility check is sensitive", {
 data(banoSEMeta) 
 bigec2 = H5S_source("http://54.174.163.77:5000")
 banoh5 = bigec2[["assays"]]
 expect_error( rr = RESTfulSummarizedExperiment(banoSEMeta[-1,], banoh5) )
})
