

load("../hdlux_all_strains_limma_fit_pvals.RData")

reordered = c(1:3,10:18,7:9,4:6)

logfc = grep("logFC",colnames(res))
logfc = res[,logfc]
logfc = logfc[,reordered]

p = grep("P.Value",colnames(res))
p = res[,p]
p = p[ , reordered ]

q = grep("adj.P.Val",colnames(res))
q = res[,q]
q = q[ , reordered ]

up = q < 0.05 & logfc > 0
dn = q < 0.05 & logfc < 0

age = c( rep( "" , 5 ) , rep("6-months",5) , rep("10-months",5) )
allele = rep( c("Q80","Q92","Q111","Q140","Q175") , 3 )

pdf("deg_hist.pdf" , width = 12 , height = 8 )
par( mfrow = c(3,6) , mar = c(2,2,2,0) , oma = c(6,7,3,0) )
for( i in 1:15 ) {
hist( 	logfc[ up[,i] , i ] ,
	breaks = -80:80 / 10 , 
	xlim = c(-2,2) ,
	ylim = c(0,1200) ,
	xlab = "" ,
	col = "green" ,
	main = "" ,
	axes = F
)
par( new = T )
hist(	logfc[ dn[,i] , i ] ,
        breaks = -80:80 / 10 ,
        xlim = c(-2,2) ,
        ylim = c(0,1200) ,
	col = "red" ,
	axes=F , main = ""
)
abline( v = 0 , lty = 2 )
if( i %in% c(2:5,7:10,12:15) ) {
axis( side = 2 , tick = F , labels = F )
}
if( i %in% c(1:10) ) {
axis( side = 1 , tick = F , labels = F )
}
if( i %in% c(1,6,11) ) {
axis( side = 2 , cex.axis = 1.8 , las = 1 )
mtext( side = 3 , line = 1 , cex = 1.5 , adj = 0 ,  age[i] )
}
if( i > 10 ) {
axis( side = 1 , cex.axis = 1.8 , las = 1 )
}
text( x = -2 , y = 1100 , cex = 2.5 , adj = 0 , allele[i] )
}
mtext( outer = T , side = 1 , line = 4 , cex = 2.5 , "logFC" )
mtext( outer = T , side = 2 , line = 4 , cex = 2.5 , "Number of Genes" )
dev.off()

ymax=max(-log10(p[,1:12]))
xmin=-1.5
xmax=1.5

p.thresh = rep(NA,18)
for( i in 1:18 ) {
  sig = which( q[,i] < 0.05 )
  p.thresh[i] = max( p[sig,i] )
}




pdf("volcano.pdf",width=7,height=10)
par( 	mfrow=c(6,3) ,
	mar = c(2,2,1,1) , cex.axis = 1 ,
	oma = c(6,6,2,0) , bty="l" )
for( i in 1:18 ) {
  smoothScatter( x = logfc[,i] ,
	y = -log10(p[,i]) ,
	xlab="", las = 1 ,
	ylab="", main="" ,
	colramp = colorRampPalette(
	  c("white","darkblue","blue","cyan","yellow","red","maroon")) ,
	xlim=c(xmin,xmax) , ylim=c(0,ymax) ,
	nrpoints = max(25,length(which(q[,i] < 0.05))) ,
	xaxt = "n" , yaxt = "n"
	)
  axis( side = 1 , at = c(-1,0,1) )
  axis( side = 2 , at = c(0,10,20) )
  if( i %in% 1:3 ) {
  	mtext( side=3 , cex = 1 , las = 0 ,
	line = 0 , 
	c("4-9 weeks" , "10-15 weeks" , "16-20 weeks")[i] )
  }
  if( i %in% c(1,4,7,10,13,16) ) {
	mtext( side = 2 , cex = 1 , line=3 , las = 0 ,
	c(rep("B6.Q111",3),
	rep("CD-1.Q111",3) ,
	rep("FVB.Q111",3) , 
	rep("129.Q111",3) ,
        rep("B6.Q92",3),
        rep("B6.Q50",3))[i] )
  }
  abline( h = -log10(p.thresh[i]) , lty = 2 )
}
mtext( side = 1 , "log2(Fold Change)" , cex= 2 , outer=T , line = 4 )
mtext( side = 2 , "-log10(p-value)" , cex=2 , outer=T , line = 4 ) 
dev.off()


















