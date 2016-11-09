## source('lib/util_agiparam.R');
source('lib/util_file.R');
source('lib/util_fit.R');
source('lib/util_metadata.R');
source('lib/util_options.R');
source('lib/util_agiparam.R');
source('lib/util_plot.R');
suppressPackageStartupMessages(library(RMySQL));

## used to set argument from defaults or parent
default=function(what,default=NULL,from='opt') {
  mget.list=mget(from,envir=parent.frame(),inherits=T,ifnotfound=list(opt=NULL));
  defaults=mget.list[[from]];
  if (is.null(defaults)) return(default);
  defaults[[what]];
}
parent=function(what) get(what,envir=parent.frame(n=2))

## not in - based on example in RefMan - must exist somewhere but...
"%notin%"=function(x,table) match(x,table,nomatch=0)==0

## get name of Rscript script. adapted from optparse source
scriptname=function() {
  script=sub("--file=","", grep("--file=",commandArgs(),value=TRUE)[1]);
  if (is.na(script)) script='interactive';
  script;
}
## project is last component of working directory. used for setting database
project=function() basename(getwd());
## connect to standard database
db_connect=function(database=project()) dbConnect(MySQL(),database);

#################### alias for browser so it'll be easier to find ####################
BROWSER=browser;
BREAKPOINT=browser;

## get index of median element of vector. analogous to which.min, which.max
which.med=function(x) {mid=ceiling(length(x)/2); which(x==sort(x)[mid])[1]}

#################### specialized versions of apply ####################
## apply f to sub-matrices of e defined by groups
gapply=function(groups,e,f) {
  t(sapply(groups,function(g) {
    eg=e[g,];
    if (length(dim(eg))>=2) apply(eg,2,function(column) f(column))
    else sapply(eg,function(column) f(column));
  }))
}
## apply f to columns of e. simple wrapper for apply
capply=function(e,f) if (!is.null(dim(e))) apply(e,2,function(column) f(column)) else f(e)
## apply f to rows of e. simple wrapper for apply
rapply=function(e,f) if (!is.null(dim(e))) apply(e,1,function(row) f(row)) else f(e)

#################### order and sample expression matrices ####################
## order expression matrix by metadata
## if order.by is NULL, assume metadata already in desired order
order_e=
  function(e=get('e',envir=parent.frame()),
           m=get('metadata',envir=parent.frame()),
           order.by=NULL) {
    if (!is.null(order.by)) m=m[order(m[[order.by]]),];
    if (length(dim(e))>=2) e=e[,rownames(m)]
      else e=e[rownames(m)];
    e;
  }
## generate random subset of rows
sample_rows=function(e,n) {
  if (length(dim(e))>=2) {
    m=dim(e)[1];
    s=sample.int(m,min(m,n));
      e=e[s,];
  } else {
    m=length(e);
    s=sample.int(m,min(m,n));
    e=e[s];
  }
  e;
}
## convert slide numbers to array numbers -- trivial but I alwyas get it wrong..
## only makes sense for 'hyb effect' order
slide.to.array.first=function(slide) 8*(slide-1)+1;
slide.to.array.last=function(slide) 8*(slide);
slide.to.array.range=function(slide.range)
  c(slide.to.array.first(slide.range[1]),slide.to.array.last(slide.range[2]));
array.to.slide=function(array) floor((array-1)/8)+1;
