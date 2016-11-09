#!/usr/bin/env Rscript
## reduce emats by aggregating replicated probes and [TODO: optionally genes]
## TODO: expr_reduce does same thing for exprs
source('lib/util.R');
## used to convert 'method' option to internal form and internal form to filename form
opt_method=list(
  mean='mean',median='median',min='min',max='max');
method_fn=list(mean='mn',median='md');

do_options=function() {
  get_options(ematdir=T,pattern.yes='^expr',pattern.no='reduce',output.prefix='reduce',method='mean'
          );
}
doit=
  function(opt=get('opt',envir=parent.frame()),
           ematdir=opt$ematdir,pattern.yes=opt$pattern.yes,pattern.no=opt$pattern.no,
           method=opt$method,
           verbose=opt$verbose) {
    method=get(method);            # convert function name to function
    suffix=c('RData','txt');
    files=findfiles(ematdir,pattern.yes,pattern.no,suffix=suffix);
    for (file in files) {
      infile=file.path(ematdir,file);
      e=readmat(infile);
      if (verbose) print(paste(sep=' ','Processing',desuffix(infile,suffix=suffix)));
      probes=rownames(e);
      probes.split=split(1:length(probes),as.factor(probes));
      ## NG 13-10-10: doing apply across entire matrix is real slow. faster to just do dups
      probes.out=names(probes.split);
      counts=sapply(probes.split,length);
      ## preallocate output matrix
      e.out=matrix(nrow=length(probes.out),ncol=ncol(e));
      colnames(e.out)=colnames(e);
      rownames(e.out)=probes.out;
      e.out[probes.out,]=e[probes.out,];
      ## '-which' below is thanks to stackoverflow!
      probes.dup=probes.split[-(which(counts==1,arr.ind=TRUE))];
      e.out[names(probes.dup),]=gapply(probes.dup,e,method);
      writedata(e.out,file=outfile());
    }
    if (opt$write.readme) write_readme(ematdir);
  }
## construct output filename
outfile=function(opt=get('opt',envir=parent.frame()),file=get('file',envir=parent.frame())) {
  if (opt$method=='mean') base=opt$output.prefix
    else base=paste(sep='_',opt$output.prefix,method2fn());
  filename(opt$ematdir,base=base,tail=file);
}

if (!interactive()) {
  do_options()
  doit();
}

