setwd("/proj/price1/sament/hdlux")

load("hdlux_all_strains_limma_fit_pvals.RData")

symbol = res$GeneSymbol

qval = grep("adj.P.Val",colnames(res))
qval = res[ , qval ]
n = ncol(qval)

deprobes = rep(NA,n)
degenes = rep(NA,n)

for( i in 1:n ) {
  x = which( qval[ , i ] < 0.05 )
  deprobes[i] = length(x)
  y = length( setdiff( unique( symbol[ x ] ) , "" ))
  degenes[i] = y
}

decounts = data.frame( 
  group = colnames(qval) ,
  deprobes , degenes )
decounts = decounts[1:18,]


