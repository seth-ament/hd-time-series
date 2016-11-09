#!/usr/bin/env Rscript
## produce quality metrics for emat
source('lib/util.R');
do_options=function() {
  get_options(ematdir=T,hmatdir=T,qualdir=T,pattern.yes='^norm',pattern.no=T,phase=2);
}
doit=
  function(opt=get('opt',envir=parent.frame()),
           ematdir=opt$ematdir,hmatdir=opt$hmatdir,qualdir=opt$qualdir,phase=opt$phase,
           pattern.yes=opt$pattern.yes,pattern.no=opt$pattern.no,
           verbose=opt$verbose) {
    createdir(qualdir);     # create output directory if necessary
    hyb_ids=get_hyb_ids(phase=phase,derived.ok=F);
    ## get stats
    suffix=c('RData','txt');
    hmatfile=filename(hmatdir,'stats');
    stats=readdata(hmatfile);
    h=hyb_ids[hyb_ids %in% colnames(stats)];
    stats=stats[,h];
    ## compute stats score: number of 'IsInRnage' attributes set
    inrange.attrs=grep('IsInRange',rownames(stats),value=T);
    inrange=stats[inrange.attrs,h];
    mode(inrange)='numeric';
    score.stats=colSums(inrange);

    files=findfiles(ematdir,pattern.yes,pattern.no,suffix=suffix);
    for (file in files) {
      infile=file.path(ematdir,file);
      e=readmat(infile);
      if (verbose) print(paste(sep=' ','Processing',desuffix(infile,suffix=suffix)));
      h=hyb_ids[hyb_ids %in% colnames(e)];
      e=e[,h];
      e.median=apply(e,1,median);
      score.cor=cor(e.median,e,method='spearman'); # correlation of arrays vs. median
      score.cor=score.cor[1,];                     # convert to vector
      score=combine_scores(score.stats,score.cor);
      q=data.frame(hyb_id=h,score.stats=score.stats,score.cor=score.cor,score=score);

      outfile=filename(qualdir,base=desuffix(file,suffix));
      writedata(q,outfile);
    }
    if (opt$write.readme) write_readme(qualdir);
  }

########## combine stats score and correlation into single score
## simple scaled Euclidean distance
combine_scores=function(stats,cor) sqrt((stats/10)^2+cor^2)/sqrt(2)

if (!interactive()) {
  do_options()
  doit();
}

