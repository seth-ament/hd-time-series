require( doBy )
require( vioplot )

load("metadata.full.RData")
 
genotype = paste( traitData$strain , traitData$allele_nominal , sep="." )

cag = traitData$allele_actual
cag[ cag == "WT" ] = 5
cag = as.numeric( gsub("Q","",cag) )
qlength = cag + 2

df = data.frame( genotype , qlength )

summ = summaryBy( cag ~ genotype , data = df , FUN = c(mean,sd) )

pdf( "qlength_distribution_vioplot.pdf" , height = 6 , width = 8 )
par( bty = "l" )
vioplot(
  df[ df[,1] == "B6.Q50" , 2 ] , 
  df[ df[,1] == "B6.Q92" , 2 ] , 
  df[ df[,1] == "B6.Q111" , 2] , 
  df[ df[,1] == "CD1.Q111" , 2] ,
  df[ df[,1] == "FVB.Q111" , 2] , 
  df[ df[,1] == "129.Q111" , 2] ,
  names = c("B6.Q50","B6.Q92","B6.Q111","CD1.Q111","FVB.Q111","129.Q111")
)
abline( h = 111 , lty = 2 )
abline( h = 92 , lty = 2 )
abline( h = 50 , lty = 2 )
segments( x0 = 1.4 , x1 = 0.6 , y0 = 50 , y1 = 50 , lwd = 3 )
mtext( side=2 , line = 3 , "HTT Q-length" )
dev.off()








