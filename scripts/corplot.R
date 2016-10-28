# scatterplots of correlations between mouse models for top 100 genes

load("../hdlux_all_strains_limma_fit_pvals.RData")

logfc = grep("logFC",colnames(res))
logfc = res[,logfc]
logfc = logfc[ , c(1:3,10:18,4:9) ]


q = grep("adj.P.Val",colnames(res))
q = res[,q]

top100 = order( res$B6.Q111_P.Value )[1:100]

a = logfc[ top100 , 3 ]
r.top100 = rep(NA , ncol(logfc) )
p.top100 = rep( NA , ncol(logfc) )
slope = rep(NA,ncol(logfc))
for( i in 1:ncol(logfc) ) {
  b = logfc[ top100 , i ]
  fit = lm( b ~ a )
  p.top100[i] = anova( fit )[[5]][1]
  slope[i] = coef( fit )[2]
  r.top100[i] = sign(slope[i]) * sqrt( summary( fit )$r.squared )
}



x = data.frame( grp = gsub("_logFC" , "" , colnames(logfc) ) , 
	r.top100 , p.top100 , slope )

write.csv( x , row.names=F , file="top100.correlations.csv" )





pdf("corplot.top100.pdf" , height = 12 , width = 5 )
par( mfrow=c(6,3) , mar = c(1,1,2,0) , oma = c(5,5,0,0) , bty = "l" )
for( i in 1:18 ) {
  a = logfc[ top100 , 3 ]
  b = logfc[ top100 , i ]
  plot( a , b , pch = 19 ,
	xlim=c(-1.5,0.7) ,
	ylim=c(-1.5,0.7) ,
	xaxt = "n" ,yaxt = "n"   
  )
  abline( lm(b~a) , lty = 1 , col = "red" )
  if( i %in% 16:18 ) {
    axis( side = 1 )
  }
  if( i %in% c(1,4,7,10,13,16) ) {
    axis( side = 2 )

    mtext( side=3 , line=0.5 , adj=-0.1 , cex=1 , 
	c("B6.Q111","CD-1.Q111","FVB.Q111","129.Q111","B6.Q92","B6.Q50")[ (i-1)/3+1 ] )
  }
  mtext( side=3 , line=-1 , adj = 0.1 , c("early","middle","late")[i] )
}
mtext( side=1 , outer=T , line=3 , "log(Fold Change) in 16-21-week-old B6.Q111 mice")
mtext( side=2 , outer=T , line=3 , "log(Fold Change) in Other Condition" )
dev.off()






