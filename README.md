# restfulSE

This R package includes proof-of-concept code illustrating several approaches to SummarizedExperiment design with
assays stored out-of-memory.

## HDF5 server backed SummarizedExperiment

[HDF Server](https://support.hdfgroup.org/projects/hdfserver/) "extends the HDF5 data model to efficiently store large data objects (e.g. up to multi-TB data arrays) and access them over the web using a RESTful API."  In this `restfulSE` package,
several data structures are introduced 

- to model the server data architecture and 
- to perform targeted extraction of numerical data from HDF5 arrays stored on the server. 

We maintain, thanks to a grant from the National Cancer Institute,
the server [http://54.174.163.77:5000/](http://54.174.163.77:5000/)
