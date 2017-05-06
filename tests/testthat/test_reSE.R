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
 expect_true(is(dataset(bigec2, "100k"), "H5S_dataset"))
 expect_true(is(bigec2[["100k"]], "H5S_dataset"))
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
 expect_true(length(tt)==6)
 expect_true(length(grep("host", tt))==4)
})

context("RESTful SE")

test_that("RESTfulSummarizedExperiment infrastructure works against server", {
 bigec2 = H5S_source("http://54.174.163.77:5000")
 data(tenx_100k_sorted)
 n100k = bigec2[["neurons100k"]]
 rr = RESTfulSummarizedExperiment( tenx_100k_sorted, n100k )
 expect_true(validObject(rr))
 expect_true(all(dim(rr)==c(26755, 100000)))
 rr2 = rr[,1:10]
 arr2 = assay(rr2)
 sumstr = structure(c(3667, 1844, 4329, 2928, 7787, 10084, 2151, 3318, 
2851, 4158), .Names = c("AAACCTGAGATAGGAG-1", "AAACCTGAGCGGCTTC-1", 
"AAACCTGAGGAATCGC-1", "AAACCTGAGGACACCA-1", "AAACCTGAGGCCCGTT-1", 
"AAACCTGAGTCCGGTC-1", "AAACCTGCAACACGCC-1", "AAACCTGCACAGCGTC-1", 
"AAACCTGCAGCCACCA-1", "AAACCTGCAGGATTGG-1"))
 expect_true(identical(apply(arr2,2,sum), sumstr))
})
 
