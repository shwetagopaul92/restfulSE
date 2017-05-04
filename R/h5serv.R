#' @import S4Vectors
#' @import SummarizedExperiment
#' @import httr
#' @import rjson
dataref = function(serverURL, host) {
paste0(serverURL, "/?host=", host)
}

setClass("H5ShContent", representation(
# "host" content
   serverURL="character",  # includes port
   host="character",
   hrefs="character",
   lastModified="character", created="character",
   root = "character"))

setClass("H5SDatasets", representation(dsuuid="character", attrs="list"), contains="H5ShContent")

setMethod("show", "H5ShContent", function(object) {
 self = object@hrefs["self"]
 cat("H5serv content for 'host'", object@host, "\n")
 cat(" at server", object@serverURL, "\n")
# cat(" use datasets() ...\n")
})
setGeneric("dim_internal", function(object, ...) standardGeneric("dim_internal"))
setMethod("dim_internal", "H5SDatasets", function(object) {
  object@attrs$shape$dims
})
#setMethod("show", "H5SDatasets", function(object) {
# cat("H5SDatasets, (1 of ", length(object@dsuuid), "): ", object@dsuuid[1], "\n", sep="")
# cat("derived from :\n")
# callNextMethod()
#})

#' extract dataset UUID and attributes from H5ShContent
#' @param object instance of H5ShContent
#' @param \dots not used yet
setGeneric("datasets", function(object, ...) standardGeneric("datasets"))
setMethod("datasets", "H5ShContent", function(object,...) {
 target = paste0(.serverURL(object), "/datasets?host=", .host(object))
 lis = fromJSON(readBin(GET(target)$content, w="character"))
 ds = new("H5SDatasets", object, dsuuid=lis$datasets)
 atts = dsAttrs(ds)
 ds@attrs=atts@attrs
 ds
})

# should write (unexported?) methods here
.serverURL = function(x) x@serverURL
.host = function(x) x@host
.root = function(x) x@root
.dsuuid = function(x, which=1) x@dsuuid[which]
setGeneric("dsAttrs", function(x, ...) standardGeneric("dsAttrs"))
setClass("H5DatasetAttributes", representation(attrs="list"), contains="H5SDatasets")
setMethod("show", "H5DatasetAttributes", function(object) {
  cat(class(object), "instance\n")
  callNextMethod()
})
setMethod("dsAttrs", "H5SDatasets", function(x, which=1, ...) {
  target = paste0(.serverURL(x), "/datasets/", .dsuuid(x, which), 
      "/attr2?host=", .host(x))
  acon = GET(target)
  ans = fromJSON(readBin(acon$content, what="character"))
  new("H5DatasetAttributes", attrs=ans, x)
})
  

#' acquire a list of REST content for an HDF5 server query, wrap in H5ShContent ('host' content)
#' @param serverURL a character string giving the URL including port number on which HDF5 server is listening
#' @param host character
#' @examples
#' # this is data served by HDF group
#' tall = getH5ShContent( serverURL = "https://data.hdfgroup.org:7258", host="tall.data.hdfgroup.org")
#' tall
#' tds = datasets(tall)
#' tds
#' tds["0:3", "0:3"]
#' tds[c(1,3),]
#' tds[,c(1,3)]
#' tds[c(1,3),c(1,3)]
getH5ShContent = function(serverURL, host) {
 stopifnot(length(serverURL)==1, length(host)==1)
 stopifnot(is(c(serverURL, host), "character"))
 basic = GET( dataref( serverURL, host ) )
 ans = fromJSON( readBin( basic$content, what = "character" ) )
 h = ans$hrefs
 nms = sapply(h, function(x) x[[2]])
 hr = sapply(h, function(x) x[[1]])
 names(hr) = nms
 new("H5ShContent", serverURL=serverURL, host=host, 
         hrefs=hr, lastModified=ans$lastModified,
         created=ans$created, root=ans$root)
}

# basic design:
# H5ShContent -- top level reference to file
# H5SDatasets -- collection of datasets within a file
#  eventually -- H5SDataset -- a specific dataset
#  Apr 28 2017 -- just assume a single dataset 
# 

setMethod("[", c("H5SDatasets", "character", "character"), function (x, i, j, ..., drop = TRUE) 
 {
 target = paste0(.serverURL(x), "/datasets/", .dsuuid(x), "/value?host=", .host(x), "&select=[", i, ",", j, "]")
 if (isTRUE(options()$verboseRSE)) message(target)
 val = GET(target)
 ans = fromJSON( readBin( val$content, what="character" ) )
 if (is.list(ans$value)) return(do.call(rbind, ans$value))
 ans$value
 })

.isH5onedim = function(x) {
 stopifnot(is(x, "H5SDatasets"))
 dims = dsAttrs(x)@attrs$shape$dims
 length(dims)==1
}

#' @import SummarizedExperiment
setMethod("[", c("H5SDatasets", "character", "missing"), function (x, i, j, ..., drop = TRUE) 
 {
 dims = dsAttrs(x)@attrs$shape$dims
# pure vector selection
 if (length(dims)==1) {
       target = paste0(.serverURL(x), "/datasets/", .dsuuid(x), 
        "/value?host=", .host(x), "&select=[", i, "]")
       val = GET(target)
       ans = fromJSON(readBin(val$content, what = "character"))
       return(ans$value)
       }
 colind = paste0("0:", as.integer(dims[2]))
 target = paste0(.serverURL(x), "/datasets/", .dsuuid(x), "/value?host=", .host(x), 
       "&select=[", i, ",", colind, "]")
 val = GET(target)
 ans = fromJSON( readBin( val$content, what="character" ) )
 (do.call(rbind, ans$value))
 })

setMethod("[", c("H5SDatasets", "missing", "character"), function (x, i, j, ..., drop = TRUE) 
 {
 stopifnot(length(j)==1)
 dims = dsAttrs(x)@attrs$shape$dims
 rowind = paste0("0:", as.integer(dims[1]))
 target = paste0(.serverURL(x), "/datasets/", .dsuuid(x), "/value?host=", .host(x), 
       "&select=[", rowind, ",", j, "]")
 val = GET(target)
 ans = fromJSON( readBin( val$content, what="character" ) )
# there may be a bug, or semantics of retrieving a column are
# different from those of retrieving a row.  a row retrieval comes
# back as an R list, a column retrieval comes back as a vector
# we need a list 
 if (!is.list(ans$value)) ans$value = list(ans$value)
 (do.call(rbind, ans$value))
 })

#trivial tuple generator
tupleGen = function(x) {
 y = x-1
 paste(y,x,1,sep=":")
}

setMethod("[", c("H5SDatasets", "numeric", "missing"), function (x, i, j, ..., drop = TRUE) {
  rowtups = tupleGen(i)
  do.call(rbind, lapply(rowtups, function(r) x[r, ]))
})

setMethod("[", c("H5SDatasets", "missing", "numeric"), function (x, i, j, ..., drop = TRUE) {
  coltups = tupleGen(j)
# retrieval of columns has different semantics
  t(do.call(rbind, lapply(coltups, function(thecol) x[, thecol ])))
})

setMethod("[", c("H5SDatasets", "numeric", "numeric"), function (x, i, j, ..., drop = FALSE) {
#
# some optimization could be done here to retrieve by column
# or row first prior to filtering the other index in R
#
# for now just retrieve rows and then filter in R
#
  rowtups = tupleGen(i)
#  following is a bad design when there are many columns
#  retrieval is very costly
# rowsel = do.call(rbind, lapply(rowtups, function(r) x[r, ]))
#  instead limit the rows retrieved, to the maximum column index
  maxcolind = paste0("0:", max(j)+1, ":1")
  rowsel = do.call(rbind, lapply(rowtups, function(r) x[r, maxcolind ]))
  rowsel[,j,drop=drop]
})

setMethod("dim", "H5SDatasets", function(x) x@attrs$shape$dim)
setMethod("show", "H5SDatasets", function(object) {
  dims = dim(object)
  cat("An HDF5-served dataset with internal dimensions", dims[1], "x", dims[2], "\n")
  okr = min(5, dims[2])
  okc = min(5, dims[1])
  cat("Northwest", okr, "x", okc, "submatrix:\n")
  print(t(object[1:okr, 1:okc]))  # consider transposing upstream when loaded into server
})

setClass("H5S_source", representation(
# "host" content
   serverURL="character", dsnames="List"))
setMethod("show", "H5S_source", function(object) {
cat("HDF5 server domain: ", object@serverURL, "\n")
cat(" There are", nrow(grps <- groups(object)), "groups.\n")
cat(" Use groups(), links(), ..., to probe and access metadata.\n")
cat(" Use dsnames() to get a list of datasets within groups.\n")
})

fixtarget = function(x) sub(".*host=(.*).hdfgroup.org", "\\1", x)

#' construct H5S_source
#' @param serverURL a URL for a port for HDF5Server
#' @param \dots not used
#' @export
H5S_source = function(serverURL, ...) {
  tmp = new("H5S_source", serverURL=serverURL, dsnames=List())
  grps = groups(tmp)
  thel = List(targs=lapply( 1:nrow(grps), 
        function(x) fixtarget(restfulSE:::hosts(links(tmp,x)))))
  tmp@dsnames = thel
  tmp
}
#' list datasetnames available in an H5S_source
#' @param src H5S_source instance
#' @export
dsnames = function(src) {
  src@dsnames
}


setClass("H5S_linkset", representation(links="list", group="character",
   source="H5S_source"))
setMethod("show", "H5S_linkset", function(object) {
cat("HDF5 server link set for group", object@group, "\n")
cat(" There are", length(object@links$links), "links.\n")
cat(" Use targets([linkset]) to extract target URLs.\n")
})

setGeneric("groups", function(object, index, ...) standardGeneric("groups"))
#' @exportMethod groups
#' @export groups
setMethod("groups", c("H5S_source", "missing"), function(object, index, ...) {
 target = paste0(.serverURL(object), "/groups")
 ans = fromJSON(readBin(GET(target)$content, w="character"))
 # find root group
 hh = t(sapply(ans$hrefs, force))
 rootgroup = sub(".*groups.", "", hh[hh[,2]=="root", 1])
 grps = c(rootgroup, unlist(ans$groups))
 nl = sapply(1:length(grps), function(x) {
   gname = grps[x]
   target = paste0(.serverURL(object), "/groups/", gname, "/links" )
   ans = fromJSON(readBin(GET(target)$content, w="character"))
   length(ans$links)
   })
 DataFrame(groups=grps, nlinks=nl)
})
setMethod("groups", c("H5S_source", "numeric"), function(object, index, ...) {
 groups(object)[index,,drop=FALSE]
})
setGeneric("links", function(object, index, ...) standardGeneric("links"))
#' @exportMethod links
#' @export links
setMethod("links", c("H5S_source", "numeric"), function(object, index, ...) {
 gname = groups(object, index)[1,1]
 target = paste0(.serverURL(object), "/groups/", gname, "/links" )
 ans = fromJSON(readBin(GET(target)$content, w="character"))
 new("H5S_linkset", links=ans, source=object, group=gname)
})
#' provide the full URLs for link members
#' @param h5linkset instance of H5S_linkset
#' @param index numeric index into link vector
#' @export
targets = function(h5linkset, index) {
 sapply(h5linkset@links$links, "[[", "target")
 }

hosts = function(h5linkset, index, cleanIP=TRUE) {
 ans = targets(h5linkset, index) #sapply(h5linkset@links$links, "[[", "target")
 ans = ans[grep("host=", ans)]
 if (cleanIP) gsub(".*host=", "host=", ans)
 else ans
 }
targetIds = function(h5linkset, index) {
 sapply(h5linkset@links$links, "[[", "id")
 }
      
setClass("H5S_dataset", representation(
  target="character", linkset="H5S_linkset"))


makeDatasetRef = function(linkset, index) {
  source = linkset@source
  thet = targets(linkset, index)
  host = sub(".*host=", "host=", thet[index])
  fromJSON(readBin(GET(paste0(source@serverURL, "/datasets?",  host=host))$content, w="character"))$datasets
}
transl = function(targ)  fromJSON(readBin(GET(targ)$content, w="character"))

setClass("H5S_datasets", representation(
  source="H5S_source", index="numeric", drop="ANY", simpleNames="character",
  shapes="list", hrefs="list", allatts="list"))
setMethod("show", "H5S_datasets", function(object) {
 cat("H5S_datasets instance:\n")
 aa = object@allatts
 alld = sapply(aa, function(x) paste(x$shape$dims, collapse=" x "))
 cr = sapply(aa, function(x) x$created)#, collapse=" x "))
 ba = sapply(aa, function(x) x$type$base)
 print(data.frame(dims=alld, created=cr, type.base=ba))
cat("---\n")
cat("Use [[ [simplename] ]] to acquire reference amenable to [i,j] subsetting.\n")
})
setClass("RESTfulH5", representation(key="character",
 dataURL="character", dataset="H5S_datasets"))
#' @exportMethod show
setMethod("show", "RESTfulH5", function(object) {
 dims = object@dataset@allatts[[ object@key ]]$shape$dims
 cat(paste(object@key, "(", paste(dims, collapse=" x "), ")"), "\n")
})

#' @exportMethod [
setMethod("[", c("RESTfulH5", "character", "character"), function(x, i, j, ..., drop=FALSE) {
#
# bracket selection passed directly to HDF5 server ... row-major
#
 uu = x@dataURL
 dims = x@dataset@allatts[[ x@key ]]$shape$dims
 ind1lims = as.numeric(strsplit(i, ":")[[1]])
 if (ind1lims[1] < 0) stop("negative starting index not allowed in i")
 if (ind1lims[2] > dims[1]) stop("i exceeds boundary for first index")
 ind2lims = as.numeric(strsplit(j, ":")[[1]])
 if (ind2lims[1] < 0) stop("negative starting index not allowed in j")
 if (ind2lims[2] > dims[2]) stop("j exceeds boundary for second index")
 uu = sub("%%SEL1%%", i, uu)
 uu = sub("%%SEL2%%", j, uu)
 val = transl(uu)$value
 if (is.list(val)) return(do.call(rbind, val))
 else return(val)
})

#' @exportMethod [[
setMethod("[[", c("H5S_datasets", "character"), function(x, i, j, ...) {
 sel = x@hrefs[[i]] # DataFrame
 sel = sel["self", 2]
 prep = sub("\\?host=", "/value?host=", sel)
 prep = paste0(prep, "&select=[%%SEL1%%,%%SEL2%%]")
 new("RESTfulH5", key=i, dataURL=prep, dataset=x)
})

#' obtain dataset metadata from a set of links
#' @param linkset instance of H5S_linkset
#' @param index numerical index of group to use
#' @param drop numerical indices of hosts to be ignored
#' @export
datasetRefs = function(linkset, index, drop=5) {
 ho = hosts(linkset, index)
 if (!is.null(drop)) ho=ho[-drop]
 source = linkset@source
 uuids = sapply(1:length(ho), function(host) {
  targ = paste0(source@serverURL, "/datasets?",  ho[host])
#  print(targ)
  fromJSON(readBin(GET(targ)$content, w="character"))$datasets
 })
 attrs = lapply(1:length(uuids), function(x) {
  targ = paste0(source@serverURL, "/datasets/", uuids[x], "?", ho[x])
  transl(targ)
  })
 HRperDS = lapply(attrs, function(x) {
    attnm = sapply(x$hrefs, "[[", 2)
    attval = sapply(x$hrefs, "[[", 1)
    ans = DataFrame(hrefName=attnm, hrefValue=attval)
    rownames(ans) = attnm
    ans
    })
 selves = sapply(HRperDS, function(x) x["self", "hrefValue"])
 simpleNames = sapply(HRperDS, function(x) sub(".*host=(.*).hdfgroup.org", "\\1", x[2,2]))
 names(attrs)=simpleNames
 names(HRperDS) = simpleNames
# list(simpleNames=simpleNames, attrs=attrs, selves=selves)
 new("H5S_datasets", source=linkset@source, index=index,
       drop=drop, simpleNames=simpleNames, shapes=list(), hrefs=HRperDS, allatts=attrs)
}

#mys = new("H5S_source", serverURL="http://54.163.220.201:5000")
#allg = groups(mys, 1)
#lin1 = links(mys,1)
