#!/usr/bin/env Rscript
## create emat by combining expr files plus annotation data
##   NOTE: all arrays must have same annotation!!
## create one emat for each combination of the 'expr" parameters, eg, ewhats
## we get annotation from data/AnnoFiles based on the arraytype
##   TODO: figure out arraytype from feparams
## save as RData structure and text files
source('lib/util.R');
do_options=function() {
  get_options(exprdir=T,ewhat=T,eparam=T,bgparam=T,offset=T,ematdir=T,output.prefix='expr',
          annodir=T,annodata=T,arraytype=T);
  }
doit=
  function(opt=get('opt',envir=parent.frame()),
           ewhats=opt$ewhat,eparams=opt$eparam,bgparams=opt$bgparam,offsets=opt$offset,
           annodata=opt$annodata,exprdir=opt$exprdir,ematdir=opt$ematdir,
           verbose=opt$verbose) {
    ## create output directory if necessary and possible
    createdir(ematdir);
    hyb_ids=get_hyb_ids(derived.ok=F);
    genes=readdata(annodata);
    ## which emats do we need?
    cases=data.frame(stringsAsFactors=F);
    if ('proc' %in% ewhats) {
      cases=rbind(cases,
        expand.grid(ewhat='proc',eparam=NA,bgparam=NA,offset=offsets,stringsAsFactors=F));
    }
    if ('raw' %in% ewhats) {
      cases=rbind(cases,
        expand.grid(ewhat='raw',eparam=eparams,bgparam=NA,offset=NA,stringsAsFactors=F));
    }
    if ('bgc' %in% ewhats) {
      cases=rbind(cases,
        expand.grid(ewhat='bgc',eparam=eparams,bgparam=bgparams,offset=offsets,stringsAsFactors=F));
    }
    ## do it
    for (i in 1:nrow(cases)) {
      case=cases[i,];
      do_one(ewhat=case$ewhat,eparam=case$eparam,bgparam=case$bgparam,offset=case$offset);
    }
    if (opt$write.readme) write_readme(ematdir);
  }
do_one=
  function(ewhat,eparam=NA,bgparam=NA,offset=NA,
           hyb_ids=get('hyb_ids',envir=parent.frame()),
           genes=get('genes',envir=parent.frame()),
           exprdir=get('exprdir',envir=parent.frame()),
           ematdir=get('ematdir',envir=parent.frame()),
           verbose=get('verbose',envir=parent.frame())) {
    e.m=matrix(nrow=dim(genes)[1],ncol=length(hyb_ids));
    colnames(e.m)=hyb_ids;
    rownames(e.m)=genes$ProbeName;

    for (hyb_id in hyb_ids) {
      exprfile=filename(exprdir,emxbase(hyb_id,ewhat,eparam,bgparam,offset));
      e=readdata(exprfile);
      e.m[,hyb_id]=e[[hyb_id]];
    }
    writedata(e.m,file=outfile(base=emxbase(NULL,ewhat,eparam,bgparam,offset)));
  }
## construct output filename
outfile=function(opt=get('opt',envir=parent.frame()),base)
  filename(opt$ematdir,base=paste(sep='_',opt$output.prefix,base));

if (!interactive()) {
  do_options()
  doit();
}

