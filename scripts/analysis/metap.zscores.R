
load("hdlux_all_strains_limma_fit_pvals.RData")

pval = grep("P.Value",colnames(res))
pval = res[ , pval[1:18] ]

fc = grep("logFC",colnames(res))
fc = res[ , fc[1:18] ]


# 1. calculate z-scores associated with the limma p-values (and fold changes)
z = matrix( NA , ncol=ncol(pval) , nrow=nrow(pval) )

z[ fc < 0 ] = qnorm( pval[ fc < 0 ] / 2 )
z[ fc > 0 ] = qnorm( (1- pval[ fc > 0 ] / 2 ) )

if( any( z == Inf ) ) {
 z[ z == Inf ] = max( z[ -which( z == Inf) ] )
}
if( any( z == -Inf ) ) {
 z[ z == -Inf ] = max( z[ -which(z==-Inf) ] )
}

colnames(z) = gsub("_P.Value","",colnames(pval))
zscores = data.frame( res[,1:2] , z )
save( zscores , file="hdlux_all_strains_limma_fit_zscores.RData" )

# 2. calculate the meta p-values by the sum of Z-scores (Stouffer's) method.
metaz = rep( NA , nrow(pval) )
metap = rep( NA , nrow(pval) )
for( i in 1:nrow(z) ) {
  z.ia = z[i,]
  z.i = sum( z.ia ) / sqrt(length(z.ia))
  metap[i] = 2*pnorm(-abs(z.i))
  metaz[i] = z.i
}
metaq = p.adjust( metap , method = "BH" )
meta = data.frame( probe = res[,1] , gene = res[,2] , metaz , metap , metaq )

meta[ order( meta$metap )[1:100] , ]

write.csv( meta , "metap.all_genes.z-score_method.csv" )


