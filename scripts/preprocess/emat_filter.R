#!/usr/bin/env Rscript
## remove bad arrays from emat
source('lib/util.R');
do_options=function() {
  get_options(ematdir=T,qualdir=T,pattern.yes='^norm',pattern.no='filter',output.prefix='filter',
          make_option('--score.cutoff',default=0.8,
                      help='Cutoff for score [default %default]'),
          make_option('--replicas',default=F,action='store_true',
                      help='Include replicas in the filtering [default %default]'),
          make_option('--best.replicas',default=T,action='store_true',
                      help='Keep replica with best quality score [default %default]'));
}
doit=
  function(opt=get('opt',envir=parent.frame()),
           ematdir=opt$ematdir,qualdir=opt$qualdir,
           pattern.yes=opt$pattern.yes,pattern.no=opt$pattern.no,
           score.cutoff=opt$score.cutoff,
           filter.replicas=opt$replicas,pick.best.replica=opt$best.replicas,
           verbose=opt$verbose) {
    metadata=get_metadata();
    suffix=c('RData','txt');
    files=findfiles(ematdir,pattern.yes,pattern.no,suffix=suffix);
    for (file in files) {
      ematfile=file.path(ematdir,file);
      qualfile=emat2qualfile(ematfile);
      e=readmat(ematfile);
      q=readdata(qualfile);
      if (verbose) print(paste(sep=' ','Processing',desuffix(ematfile,suffix=suffix)));
      m=merge(metadata,q,all.x=T);
      rownames(m)=m$hyb_id; # merge trashes rownames
      m=subset(m,subset=(hyb_id %in% colnames(e)));
      ## m=metadata[colnames(e),];  
      ## rows with score NA are derived from replicas. always keep these
      if (filter.replicas) m=subset(m,subset=(is.na(score)|score>=score.cutoff))
        else m=subset(m,subset=(is.na(score)|score>=score.cutoff|!is.na(replica_number)));
      if (pick.best.replica) {
        m.replicas=subset(m,subset=(!is.na(replica_number)));
        g.replicas=split(m.replicas,m.replicas$sample_id);
        ## the do.call(rbind,..) construct below is from R Cookbook Recipe 5.19
        best.replicas=do.call(rbind,lapply(g.replicas,function(g) g[which.max(g$score),]))
        rownames(best.replicas)=best.replicas$hyb_id;
        m=rbind(subset(m,subset=(is.na(replica_number))),best.replicas);
      }
      e=e[,m$hyb_id];
      writedata(e,outfile());
    } 
    if (opt$write.readme) write_readme(ematdir);
  }
## construct output filename
outfile=function(opt=get('opt',envir=parent.frame()),file=get('file',envir=parent.frame()))
  filename(opt$ematdir,base=opt$output.prefix,tail=file);

if (!interactive()) {
  do_options();
  doit();
}
