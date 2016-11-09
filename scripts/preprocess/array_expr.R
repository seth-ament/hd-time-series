#!/usr/bin/env Rscript
## read Agilent array files, extract expression, and save as RData structures and text files
suppressPackageStartupMessages(library(limma));
suppressPackageStartupMessages(library(optparse));
source('lib/util.R');
do_options=function() {
    get_options(arraydir=T,pattern.yes='hyb',exprdir=T,ewhat=T,eparam=T,bgparam=T,offset=T);
}
doit=
  function(opt=get('opt',envir=parent.frame()),
           ewhats=opt$ewhat,eparams=opt$eparam,bgparams=opt$bgparam,offsets=opt$offset,
           arraydir=opt$arraydir,exprdir=opt$exprdir,pattern.yes=opt$pattern.yes,
           verbose=opt$verbose) {
    ## create output directory if necessary and possible
    createdir(exprdir);
    ## which exprs do we need to make? code from expr_emat
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
    ## make exprs for each array file. only make file if it doesn't exist
    arrayfiles=findfiles(arraydir,pattern.yes,suffix='txt');
    for (arrayfile in arrayfiles) {
      hyb_id=desuffix(arrayfile,suffix='txt');
      exprfiles=apply(cases,1,function(case)
        exprfile=filename(exprdir,
          emxbase(hyb_id,case['ewhat'],case['eparam'],case['bgparam'],case['offset'])));
      need=sapply(exprfiles,function(exprfile) !filesexist(exprfile));      
      skip=which(!need);
      if (verbose&&length(skip)>0)
        print(paste(sep='','Skipping ',paste(collapse=', ',exprfiles[skip])));
      need=which(need);
      ## read elists we need
      if (length(need)==0) next;
      if (verbose) print(paste(sep='','Reading ',filename(arraydir,arrayfile)));
      elists<<-list();
      ewhats=cases[need,'ewhat'];
      if ('proc' %in% ewhats)
        elists$proc<<-read_array(eparam='proc',bgparam=NA); # other args set automatically
      ## if 'bgc' present, defines all remaining elists we need. else try 'raw'
      params=unique(subset(cases[need,],select=c(eparam,bgparam),subset=(ewhat=='bgc')))
      if (nrow(params)==0)
        params=unique(subset(cases[need,],select=c(eparam,bgparam),subset=(ewhat=='raw')));
      if (nrow(params)>0) 
        for (i in 1:nrow(params)) {
          eparam=params[i,'eparam'];
          bgparam=params[i,'bgparam'];
          elist=read_array();       # args set automatically
          ## use for raw and bgc
          elists[[eparam]]<<-elist;                           # for raw
          if (!is.na(bgparam)) elists[[paste(sep='.',eparam,bgparam)]]<<-elist;   # for bgc
        }
      for (i in need) {
        case=cases[i,];
      do_one(ewhat=case$ewhat,eparam=case$eparam,bgparam=case$bgparam,offset=case$offset);
      }}
    if (opt$write.readme) write_readme(exprdir);
  }

do_one=
  function(ewhat,eparam=NA,bgparam=NA,offset=NA,
           hyb_id=get('hyb_id',envir=parent.frame()),
           exprdir=get('exprdir',envir=parent.frame()),
           verbose=get('verbose',envir=parent.frame())) {
    exprfile=filename(exprdir,emxbase(hyb_id,ewhat,eparam,bgparam,offset));
    offset=as.numeric(offset);
    if (ewhat=='proc') {
      elist=elists$proc;
      write_expr(elist$E+offset);
    } else if (ewhat=='raw') {
      elist=elists[[eparam]];                 # args set automatically
      write_expr(elist)
    } else {
      ## ewhat must be bgc
      elist=elists[[paste(sep='.',eparam,bgparam)]];
      elist.bgc=backgroundCorrect(elist,method="normexp",offset=offset,verbose=F);
      write_expr(elist.bgc);
    }
  }
read_array=
  function(arrayfile=get('arrayfile',envir=parent.frame()),
           arraydir=get('arraydir',envir=parent.frame()),
           eparam=get('eparam',envir=parent.frame()),
           bgparam=get('bgparam',envir=parent.frame()),
           verbose=get('verbose',envir=parent.frame())) {
    eparam=opt2param(eparam,'fg');
    bgparam=opt2param(bgparam,'bg');
    if (is.na(bgparam))
      read.maimages(files=arrayfile,path=arraydir,
                    columns=list(E=eparam),annotation=c("ControlType","ProbeName"),
                    verbose=F)
      else
        read.maimages(files=arrayfile,path=arraydir,
                      columns=list(E=eparam,Eb=bgparam),annotation=c("ControlType","ProbeName"),
                      verbose=F);
  }
write_expr=
  function(e,
           exprfile=get('exprfile',envir=parent.frame()),
           hyb_id=get('hyb_id',envir=parent.frame()),
           verbose=get('verbose',envir=parent.frame())) {
    if (mode(e)=='list') e=e$E;
    e=data.frame(e);
    colnames(e)=hyb_id;
    writedata(e,exprfile,row.names=F,verbose=verbose);
  }
## write_annodata=
##   function(genes,file='annodata',
##            exprdir=get('exprdir',envir=parent.frame()),
##            verbose=get('verbose',envir=parent.frame())) {
##     annodata=file.path(exprdir,file);
##     if (verbose) print(paste(sep=' ','Writing',annodata));
##     save(genes,file=paste(sep='.',annodata,'RData'));
##     writetxt(genes,file=paste(sep='.',annodata,'txt'));
##   }
## write_readme=
##   function(opt,file='README',
##            exprdir=get('exprdir',envir=parent.frame()),
##            verbose=get('verbose',envir=parent.frame())) {
##     lines=NULL;
##     lines=c(lines,paste(sep=' ','created',date()));
##     for (param in names(opt)) {
##       lines=c(lines,paste(sep='=',param,paste(collapse=', ',opt[[param]])));
##     }
##     file=file.path(exprdir,file);
##     if (verbose) print(paste(sep=' ','Writing',file));
##     writeLines(lines,con=file);
##   }
if (!interactive()) {
  do_options()
  doit();
}
