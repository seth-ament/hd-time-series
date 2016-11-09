#!/usr/bin/env Rscript
## normalize emats in a specified directory
suppressPackageStartupMessages(library(limma));
source('lib/util.R');
## used to convert 'method' option to internal form and internal form to filename form
opt_method=list(
  quantile='quantile',qntl='quantile',
  cyclicloess='cyclicloess',loess='cyclicloess',
  scale='scale');
method_fn=list(quantile='qntl',cyclicloess='loess',scale='scale');

do_options=function() {
  get_options(ematdir=T,pattern.yes='^reduce',pattern.no='norm',output.prefix='norm',
          method='quantile',
          make_option('--log2.transform',default=T,action='store_true',
                      help='log2 transform data after normalization [default %default]'));
}
doit=
  function(opt=get('opt',envir=parent.frame()),
           ematdir=opt$ematdir,pattern.yes=opt$pattern.yes,pattern.no=opt$pattern.no,
           method=opt$method,log2.transform=opt$log2.transform,
           verbose=opt$verbose) {
    method=match.arg(method,c('quantile','cyclicloess','scale'));
    if (method=='loess') method='cyclicloess';
    suffix=c('RData','txt');
    files=findfiles(ematdir,pattern.yes,pattern.no,suffix=suffix);
    for (file in files) {
      infile=file.path(ematdir,file);
      e=readmat(infile);
      if (verbose) print(paste(sep=' ','Processing',desuffix(infile,suffix=suffix)));
      e.out=normalizeBetweenArrays(e,method=method);
      ## note that limma does not log-transform matrices
      if (log2.transform) e.out=log2(e.out);
      writedata(e.out,file=outfile());
    }
    if (opt$write.readme) write_readme(ematdir);
  }
## construct output filename
outfile=function(opt=get('opt',envir=parent.frame()),file=get('file',envir=parent.frame())) {
  if (opt$method=='quantile') base=opt$output.prefix
    else base=paste(sep='_',opt$output.prefix,method2fn());
  filename(opt$ematdir,base=base,tail=file);
}

if (!interactive()) {
  do_options();
  doit();
}

