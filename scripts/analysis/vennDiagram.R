

load("../hdlux_all_strains_limma_fit_pvals.RData")

q = grep("adj.P.Val",colnames(res))
q = res[,q]

degs = as.matrix(q)
degs[ q < 0.05 ] = 1
degs[ q >= 0.05 ] = 0

genes = setdiff( unique( res[,2] ) , c(NA,"") )
b6.q111.middle = as.numeric( genes %in% unique( res[ q[,2] < 0.05 , 2 ] ) )
b6.q111.late = as.numeric( genes %in% unique( res[ q[,3] < 0.05 , 2 ] ) )
cd1.q111.late = as.numeric( genes %in% unique( res[ q[,12] < 0.05 , 2 ] ) )

degs = data.frame( b6.q111.middle , b6.q111.late , cd1.q111.late )
rownames(degs) = genes

library(limma)

pdf("venn.pdf")
vennDiagram( degs )
dev.off()

t = table( degs[,1] == 1 , degs[,2] == 1 )
fisher.test( t , alternative = "greater" )$p.value
# p = 1.0e-174

t = table( degs[,1] == 1 , degs[,3] == 1 )
fisher.test( t , alternative = "greater" )$p.value
# p = 2.0e-94

t = table( degs[,2] == 1 , degs[,3] == 1 )
fisher.test( t , alternative = "greater" )$p.value
# p = 1.3e-173



