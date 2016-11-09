
load("../go_enrichments.fisher.limma_p001.RData")

load("../hdlux_all_strains_limma_fit_pvals.RData")

p.gene = grep("P.Value",colnames(res))
p.gene = res[,p.gene]


q.gene = grep("adj.P.Val",colnames(res))
q.gene = res[ , q.gene ]

a = which( q.gene[,1] < 0.05 )
b = which( q.gene[,2] < 0.05 )
c = which( q.gene[,3] < 0.05 )

genes = res[,2]

sort( genes[ intersect( b , c ) ] )

p.path.dn = fisher.res$down$p
q.path.dn = apply( p.path.dn , 2 , p.adjust , method = "BH" )

colSums( q.path.dn < 0.05 )

t = table( q.path.dn[,2] < 0.05 , q.path.dn[,3] < 0.05 )
fisher.test( t , alternative = "greater" )

paths = rownames( p.path.dn )
paths[ q.path.dn[,2] < 0.05 & q.path.dn[,3] < 0.05 ]


