#' H5S_source identifies an HDF5 server and manages some metadata about contents
#' 
#' @name H5S_source
#' @rdname H5S_source-class
#' @slot serverURL character string with a URL
#' @slot dsmeta DataFrame instance with metadata about content of server
#' @exportClass H5S_source
setClass("H5S_source", representation(
# "host" content
   serverURL="character", dsmeta="DataFrame"))
setMethod("show", "H5S_source", function(object) {
cat("HDF5 server domain: ", object@serverURL, "\n")
cat(" There are", nrow(object@dsmeta), "groups.\n")
cat(" Use groups(), links(), ..., to probe and access metadata.\n")
cat(" Use dsmeta() to get information on datasets within groups.\n")
cat(" Use [[ [dsname] ]]  to get a reference suitable for [i, j] subsetting.\n")
})

#' name H5S_dataset
#' rdname H5S_dataset-class
#' @slot source instance of H5S_source instance
#' @slot simpleName character string naming dataset 
#' @slot shapes list including dimension information
#' @slot hrefs DataFrame of hrefs as defined in the API
#' @slot allatts list of all attributes
#' @slot presel string prepared for select operation in GET
#' @exportClass H5S_dataset
setClass("H5S_dataset", representation(
  source="H5S_source", simpleName="character",
  shapes="list", hrefs="DataFrame", allatts="list", presel="character"))
setMethod("show", "H5S_dataset", function(object) {
 cat("H5S_dataset instance:\n")
# aa = object@allatts
# alld = sapply(aa, function(x) paste(x$shape$dims, collapse=" x "))
# cr = sapply(aa, function(x) x$created)#, collapse=" x "))
# ba = sapply(aa, function(x) x$type$base)
 curdim = object@shapes$dims
 print(data.frame(dsname=object@simpleName, intl.dim1=curdim[1], intl.dim2=curdim[2], 
     created=object@allatts$created, type.base=object@allatts$type$base))
#cat("Use [[", object@simplename, "]] to acquire reference amenable to [i,j] subsetting.\n")
})


fixtarget = function(x) sub(".*host=(.*).hdfgroup.org", "\\1", x)

#' construct H5S_source
#' @name H5S_source
#' @rdname H5S_source-class
#' @param serverURL a URL for a port for HDF5Server
#' @param \dots not used
#' @note The dsmeta slot holds a DataFrame with a column \code{dsnames}
#' that is a list with ith element a character vector of all dsnames
#' available for the ith group.  There is no effort at present to
#' search all groups for candidate datasets.
#' @examples
#' bigec2 = H5S_source("http://54.174.163.77:5000")
#' bigec2
#' dsmeta(bigec2)[1:2,] # two groups
#' dsmeta(bigec2)[1,2][[1]] # all dataset candidates in group 1
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
#' @rdname H5S_source-class
#' @exportMethod [[
setMethod("[[", c("H5S_source", "character", "ANY"), function(x, i, j) {
  dataset(x, i)
})



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
      


makeDatasetRef = function(linkset, index) {
  source = linkset@source
  thet = targets(linkset, index)
  host = sub(".*host=", "host=", thet[index])
  fromJSON(readBin(GET(paste0(source@serverURL, "/datasets?",  host=host))$content, w="character"))$datasets
}
transl = function(targ)  fromJSON(readBin(GET(targ)$content, w="character"))

# SHOULD BE DELETED
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
# PROBABLY NOT NEEDED
setClass("RESTfulH5", representation(key="character",
 dataURL="character", dataset="H5S_datasets"))
#' @exportMethod show
setMethod("show", "RESTfulH5", function(object) {
 dims = object@dataset@allatts[[ object@key ]]$shape$dims
 cat(paste(object@key, "(", paste(dims, collapse=" x "), ")"), "\n")
})

#' @name H5S_dataset
#' @rdname H5S_source-class
#' @exportMethod [
setMethod("[", c("H5S_dataset", "character", "character"), function(x, i, j, ..., drop=FALSE) {
#
# bracket selection passed directly to HDF5 server ... row-major
#
 uu = x@presel
 dims = x@shapes$dims
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


#' @name dataset
#' @rdname H5S_source-class
#' @param h5s instance of H5S_source
#' @param tag character string identifying a dataset
#' @export
dataset = function(h5s, tag) {
 dsns = dsmeta(h5s)[,"dsnames"]
 hits = sapply(dsns, function(x) length(grep(tag, x))>0)
 if (!any(hits)) stop("tag not found in dsmeta(h5s)")
 if (sum(hits)>1) warning("tag occurs in several groups, using first")
 #grp = dsns[which(hits),,drop=FALSE]
 fulldsn = dsns[which(hits)][[1]] # unlist
 fulldsn = fulldsn[ grep(tag, fulldsn) ]
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
 attrs = transl( sub("datasets", paste0("datasets/", uuid), targ ) )
 hrnm = sapply(attrs$hrefs, "[[", 2)
 hrval = sapply(attrs$hrefs, "[[", 1)
 ans = DataFrame(hrefName=hrnm, hrefValue=hrval)
 rownames(ans) = hrnm
 self = ans["self", "hrefValue"]
 prep = sub("\\?host=", "/value?host=", self)
 prep = paste0(prep, "&select=[%%SEL1%%,%%SEL2%%]")
 list(uuid=uuid, hrefs=ans, attrs=attrs)
 new("H5S_dataset", source=h5s, simpleName=fulldsn,
    shapes=attrs$shape, hrefs=ans, allatts=attrs, presel=prep)
}


#mys = new("H5S_source", serverURL="http://54.163.220.201:5000")
#allg = groups(mys, 1)
#lin1 = links(mys,1)
