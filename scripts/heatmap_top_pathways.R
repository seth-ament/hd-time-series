# create a heatmap with -log10(p-values) of top pathways
# Seth Ament
# 4/28/2016

load("../go_enrichments.fisher.limma_p001.RData")

# down-regulated pathways

p = fisher.res$down$p
reordered = c(1:3,10:18,4:9)
p = p[,reordered]

library( metap )
metap = rep(NA,nrow(p))
for( i in 1:nrow(p) ) {
  metap[i] = sumlog( p[i,] )$p
}

top50 = order( metap )[1:50]

x = -log10(p)

my_palette = colorRampPalette( c("white","blue") )(299)
breaks = c( seq( 0 , 7 , length=299 ) ,
		max(x[top50,]) )

pdf("heatmap.top50.pathways.logp.pdf" )
heatmap(
	x[top50,] ,
	col=my_palette ,
	breaks = breaks ,
	Colv = NA ,
	scale = "none" ,
	margins = c(5,15) )
dev.off()

# up-regulated pathways

p = fisher.res$up$p
reordered = c(1:3,10:18,4:9)
p = p[,reordered]

library( metap )
metap = rep(NA,nrow(p))
for( i in 1:nrow(p) ) {
  metap[i] = sumlog( p[i,] )$p
}

top50 = order( metap )[1:50]

x = -log10(p)

my_palette = colorRampPalette( c("white","blue") )(299)
breaks = c( seq( 0 , 7 , length=299 ) ,
                max(x[top50,]) )

pdf("heatmap.top50.up.pathways.logp.pdf" )
heatmap(
        x[top50,] ,
        col=my_palette ,
        breaks = breaks ,
        Colv = NA ,
        scale = "none" ,
        margins = c(5,15) )
dev.off()



