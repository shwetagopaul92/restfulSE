
setClass("H5S_source", representation(
# "host" content
   serverURL="character", dsmeta="DataFrame"))
setMethod("show", "H5S_source", function(object) {
cat("HDF5 server domain: ", object@serverURL, "\n")
cat(" There are", nrow(grps <- groups(object)), "groups.\n")
cat(" Use groups(), links(), ..., to probe and access metadata.\n")
cat(" Use dsmeta() to get information on datasets within groups.\n")
})

fixtarget = function(x) sub(".*host=(.*).hdfgroup.org", "\\1", x)

#' construct H5S_source
#' @param serverURL a URL for a port for HDF5Server
#' @param \dots not used
#' @examples
#' bigec2 = H5S_source("http://54.174.163.77:5000")
#' bigec2
#' dsmeta(bigec2)[1:2,]
#' @export
H5S_source = function(serverURL, ...) {
  tmp = new("H5S_source", serverURL=serverURL, dsmeta=DataFrame())
  grps = groups(tmp)
  message("analyzing groups for their links...")
  thel = targs=List(targs=lapply( 1:nrow(grps), 
        function(x) fixtarget(restfulSE:::hosts(links(tmp,x)))))
  message("done")
  tmp@dsmeta = DataFrame(groupnum=1:nrow(grps), dsnames=thel, grp.uuid=grps$groups)
  tmp
}
#' list information about datasets available in an H5S_source
#' @param src H5S_source instance
#' @export
dsmeta = function(src) {
  src@dsmeta
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

dataset = function(h5s, tag) {
 dsns = dsmeta(h5s)[,"dsnames"]
 hits = sapply(dsns, function(x) length(grep(tag, x))>0)
 if (!any(hits)) stop("tag not found in dsmeta(h5s)")
 if (sum(hits)>1) warning("tag occurs in several groups, using first")
 grp = dsns[which(hits),,drop=FALSE]
 lin = links(h5s, which(hits))
 targs = targets(lin)
 targs = targs[grep(tag,targs)]
 if (length(targs)>1) {
      warning("dataset tag does not identify a single target.  found:")
      print(targs)
      stop("please supply tag to identify single target.")
      }
 targ = sub(".host", "datasets?host", targs)
 uuid = fromJSON(readBin(GET(targ)$content, w="character"))$datasets
 attrs = transl( sub("datasets", paste0("datasets/", uuid), targs ) )
 attrs
}


#mys = new("H5S_source", serverURL="http://54.163.220.201:5000")
#allg = groups(mys, 1)
#lin1 = links(mys,1)
