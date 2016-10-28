# create a heatmap with z-scores of top genes
# Seth Ament
# 4/28/2016

load("../hdlux_all_strains_limma_fit_zscores.RData")

zmat = as.matrix( zscores[,-c(1:2)] )
rownames(zmat) = zscores[,2]

reordered = c(1:3,10:18,4:9)

zmat = zmat[,reordered]

options(stringsAsFactors=F)

metaz = read.csv("../metap.all_genes.z-score_method.csv" )

top100 = order( abs(zmat[,"B6.Q111.late"] ) , decreasing = T )[1:100]

my_palette = colorRampPalette( c("blue","white","orange") )(299)
breaks = c( 	min(zmat) ,
		seq( -7 , 0 , length=149 ) ,
		seq( 0 , 7 , length=149 ) ,
		max(zmat) )

pdf("heatmap.top100.zscores.pdf" )
heatmap(
	zmat[top100,] ,
	col=my_palette ,
	breaks = breaks ,
	Colv = NA ,
	scale = "none" , cexRow = 0.4 ,
	margins = c(5,18) )
dev.off()




