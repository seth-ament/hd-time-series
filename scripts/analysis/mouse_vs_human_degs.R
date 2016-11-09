options(stringsAsFactors=F)
setwd("/proj/price1/sament/hdlux")


# load the mouse limma results
load("hdlux_all_strains_limma_fit_pvals.RData")
limma.mouse = res
rm(res)
 
fdr = grep( "adj.P.Val" , names(limma.mouse) )
fdr = limma.mouse[ , fdr ]

logfc = grep( "logFC" , names(limma.mouse) )
logfc = limma.mouse[ , logfc ]

mouse.genes = limma.mouse$GeneSymbol

# load the human limma results

load("../hodges/trn_2015-8-2/limma_pvals.RData")
limma.human = res
rm(res)

human.dn = limma.human[ limma.human$adj.P.Val < 0.05 & limma.human$logFC < 0 , "hgnc_symbol" ]
human.dn = unique( human.dn )
human.up = limma.human[ limma.human$adj.P.Val < 0.05 & limma.human$logFC > 0 , "hgnc_symbol" ]
human.up = unique( human.up )

# 1:1 orthologs based on MGI annotations (downloaded 12/16/2015)
orth = read.delim("../resources/HOM_MouseHumanSequence.rpt")
orth.m = orth[ orth$NCBI.Taxon.ID == 10090 , c(1,4) ]
orth.h = orth[ orth$NCBI.Taxon.ID == 9606 , c(1,4) ]
orth.hm = merge( orth.m , orth.h , by=1 )
dup = which( duplicated(orth.hm[,2]) ==T | duplicated(orth.hm[,3]) == T )
orth = orth.hm[-dup,]

fisher.p = matrix( 1 , nc = 2 , nr = 18 )
fisher.o = matrix( 0 , nc = 2 , nr = 18 )
fisher.or = matrix( 0 , nc = 2 , nr = 18 )

for( i in 1:2 ) {
 if( i == 1 ) {
  human.degs = human.dn
 }
 if( i == 2 ) {
  human.degs = human.up
 }
 for( j in 1:18 ) {
  if( i == 1 ) {
   mouse.degs = mouse.genes[ fdr[,j] < 0.05 & logfc[,j] < 0 ]
   mouse.degs = setdiff( unique( mouse.degs ) , "" )
  }
  if( i == 2 ) {
   mouse.degs = mouse.genes[ fdr[,j] < 0.05 & logfc[,j] > 0 ]
   mouse.degs = setdiff( unique( mouse.degs ) , "" )
  }
  human.deg.orth = which( orth[,3] %in% human.degs )
  mouse.deg.orth = which( orth[,2] %in% mouse.degs )
  t = table( 1:nrow(orth) %in% human.deg.orth , 1:nrow(orth) %in% mouse.deg.orth )
  if( any( dim(t) != 2 ) ) next
  test = fisher.test( t , alternative = "greater" )
  fisher.p[j,i] = test$p.value
  fisher.o[j,i] = t[2,2]
  fisher.or[j,i] = test$estimate
 }
}

colnames(fisher.p) = c("dn","up")
rownames(fisher.p) = gsub( "_adj.P.Val" , "" , colnames(fdr) )
fisher.p

fdr05.dn = colSums( fdr < 0.05 & logfc < 0 )
fdr05.up = colSums( fdr < 0.05 & logfc > 0 )

hmax = max( fdr05.up )
hmin = max( fdr05.dn )

library( gplots )
pdf( "fdr05.counts.barplot.pdf" )

barplot2( 	height = matrix( fdr05.up , nrow=3 ) ,
		beside=T ,
		ylim=c(-hmin,hmax) ,
		col="black"	)
par( new=T )

barplot2(       height = matrix( -fdr05.dn , nrow=3 ) ,
                beside=T ,
                ylim=c(-hmin,hmax) ,
                col="black"     )
par( new=T )

barplot2(       height = matrix( fisher.o[,2] , nrow=3 ) ,
                beside=T ,
                ylim=c(-hmin,hmax) ,
                col="grey"     )
par( new=T )

barplot2(       height = matrix( -fisher.o[,1] , nrow=3 ) ,
                beside=T ,
                ylim=c(-hmin,hmax) ,
                col="grey"     )

dev.off()





















