setwd("/proj/price1/sament/hdlux")
print(load("metadata.full.RData"))

m = traitData

group = paste( m$strain , m$allele_nominal , sep="." )

mtx = as.matrix( table( group , m$week ) )

mtx2 = matrix( 0 , nrow = 13 , ncol = 17 )
mtx2[ c( 1:4 , 6:7 , 9:10 , 12:13 ) , ] = mtx[ c(6,4,5,3,9,7,11,10,2,1) , ]
colnames(mtx2) = paste( "wk", colnames(mtx) , sep ="" )
rownames(mtx2)[c( 1:4 , 6:7 , 9:10 , 12:13 )] = rownames(mtx)[ c(6,4,5,3,9,7,11,10,2,1)]

my_palette = colorRampPalette( c("white" , "blue") )(11)

pdf("study_design_bubbleplot.pdf" , width = 8 , height = 4 )
par( bty = "l" , mar = c(4,8,3,1) )
plot( x = c(4,20) , y = c(0,14) , type = "n" , xlab = "Age (Weeks)" , 
   ylab = "" , yaxt = "n" )
points( x = sort( rep( 4:20 , 13 ) ), y = rep(13:1,17) , 
   cex = as.vector( mtx2 )/3 , pch=21 , 
   bg = my_palette[ as.vector(mtx2)+1 ] )

points( xpd = NA , x = 4:8 , y = rep(16,5), 
   cex = 1:5/1.5 , pch=21 , bg = my_palette[ c(3,5,7,9,11) ] )
text( xpd = NA , x = 4:8 , y = rep(14.5,5) , 1:5*2 )
axis( side = 2 , las = 1 , at = 13:1 , labels = rownames(mtx2) )
dev.off()





