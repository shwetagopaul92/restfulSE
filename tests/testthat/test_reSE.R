library(restfulSE)
library(rhdf5client)
library(restfulSEData)

context("connection")

test_that("H5S_source completes", {
 bigec2 = H5S_source("http://54.174.163.77:5000")
 expect_true(is(bigec2, "H5S_source"))
})

context("content wrapper structure") 

test_that("H5S_source processes", {
 bigec2 = H5S_source("http://54.174.163.77:5000")
 expect_true(all(dim(rhdf5client::groups(bigec2))==c(10,2))) 
 expect_true(is(rhdf5client::links(bigec2,1), "H5S_linkset"))
 expect_true(is(dataset(bigec2, "tenx_100k"), "H5S_dataset"))
 expect_true(is(bigec2[["tenx_100k"]], "H5S_dataset"))
 expect_true(is(dsmeta(bigec2), "DataFrame"))
 expect_true(all(dim(dsmeta(bigec2))==c(10,3))) 
})

context("indexing infrastructure")

test_that("sproc/isplit work", {
 expect_true(length(isplit(c(1,2,3,4,5,10,15,20,30:40)))==3)
 ii = isplit(c(1,2,3,4,5,10,15,20,30:40))
 ss = structure(c("0:5:1", "9:20:5", "29:40:1"), .Names = c("1", 
"2", "3"))
 expect_true(identical(ss, unlist(sproc(ii))))
 ii = isplit(c(1:10, seq(50,25,-5), seq(80,100,2)))
 ss = structure(c("0:10:1", "50:24:-5", "79:100:2"), 
                .Names = c("1", "2", "3"))
 expect_true(identical(ss, unlist(sproc(ii))))
})

context("targets generation")

test_that("targets method works", {
 bigec2 = H5S_source("http://54.174.163.77:5000")
 tt = targets(rhdf5client::links(bigec2, 1))
 expect_true(length(tt)>=9)  # july 5 2017
 expect_true(length(grep("host", tt))>=7) # added tenx_full dataset
})

context("RESTful SE")

test_that("RESTfulSummarizedExperiment infrastructure works against server", {
 bigec2 = H5S_source("http://54.174.163.77:5000")
 #data("st100k", package="restfulSEData")
 ehub = ExperimentHub::ExperimentHub()
 myfiles <- AnnotationHub::query(ehub , "restfulSEData")
 myfiles[["EH552"]] -> st100k
 n100k_served = bigec2[["tenx_100k_sorted"]]  
 rr = RESTfulSummarizedExperiment( st100k, n100k_served )
 expect_true(is(rr, "RangedSummarizedExperiment"))
 #data("banoSEMeta", package="restfulSEData")
 myfiles[["EH551"]] -> banoSEMeta
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
# library(yriMulti)
# data(banovichSE)
# bansel = assay(banovichSE[c(1,3,5,200000,300000),c(1,5,7)])
 bansel = structure(c(0.473396289, -0.222581829, -0.150630834, 1.754529575, 
-0.604892835, -0.903498324, -0.027544132, 0.596699897, 0.317356179, 
0.393075519, 1.377916164, -1.473102055, -0.268610844, 0.63228941, 
-0.286316831), .Dim = c(5L, 3L), .Dimnames = list(c("cg00000029", 
"cg00000236", "cg00000363", "cg16010467", "cg25249241"), c("NA18498", 
"NA18516", "NA18519")))
 library(restfulSE)
 #data("banoSEMeta", package="restfulSEData") 
 ehub = ExperimentHub::ExperimentHub()
 myfiles <- AnnotationHub::query(ehub , "restfulSEData")
 myfiles[["EH551"]] -> banoSEMeta
 bigec2 = H5S_source("http://54.174.163.77:5000")
 banh = bigec2[["assays"]]
 rr = RESTfulSummarizedExperiment( banoSEMeta, banh )
 options(scipen = 999)
 rrsel = rr[c(1,3,5,200000,300000),c(1,5,7)]
 rrass = assay(rrsel)
 expect_true( max(abs(rrass-bansel))<1e-6 )
 rrsel0 = rr[c(1,3,5,200000,300000),c(1)]
 rrass0 = assay(rrsel0)
 expect_true(ncol(rrass0)==1)
})

test_that("dim compatibility check is sensitive", {
 #data("banoSEMeta", package="restfulSEData") 
 ehub = ExperimentHub::ExperimentHub()
 myfiles <- AnnotationHub::query(ehub , "restfulSEData")
 myfiles[["EH551"]] -> banoSEMeta
 bigec2 = H5S_source("http://54.174.163.77:5000")
 banoh5 = bigec2[["assays"]]
 expect_error(rr==RESTfulSummarizedExperiment(banoSEMeta[-1,], banoh5) )
})
