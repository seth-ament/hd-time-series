# compute enrichments of DEGs for GO terms
options(stringsAsFactors=F)


# load the limma results 
load( "hdlux_all_strains_limma_fit_pvals.RData" )
pval = grep( "P.Value" , names(res) )
pval = res[ , pval ]
logfc = grep("logFC" , names(res) )
logfc = res[ , logfc ]
fdr = grep( "adj.P.Val" , names(res) )
fdr = res[ , fdr ]
all.genes = res$GeneSymbol
ncond = ncol( pval )

# load GO term annotations
load( "../resources/goterms_mmu_2015-11-4.RData")
nset = length(genesets)
size = rep( NA , nset ) 
for( i in 1:nset) {
 size[i] = length(genesets[[i]])
}

sets = genesets[ size >= 10 & size <= 500 ]
nset = length(sets)

fisher.res = list()

for( k in 1:3 ) {

cat( "Working on" , k , "\n" )

p = matrix( 1 , nr = nset , nc = ncond )
o = matrix( 0 , nr = nset , nc = ncond )
or = matrix( 0 , nr = nset , nc = ncond )

colnames(p) = colnames(o) = colnames(or) = gsub( "\\_(.*)" , "" , colnames(pval) )
rownames(p) = rownames(o) = rownames(or) = names(sets)

for( i in 1:ncond ) {
 if( k == 1 ) {
  degs = all.genes[ fdr[,i] < 0.05 ]
  degs = setdiff( unique(degs) , "" )
 }
 if( k == 2 ) {
  degs = all.genes[ fdr[,i] < 0.05 & logfc[,i] < 0 ]
  degs = setdiff( unique(degs) , "" )
 }
 if( k == 3 ) {
  degs = all.genes[ fdr[,i] < 0.05 & logfc[,i] > 0]
  degs = setdiff( unique(degs) , "" )
 }

 for( j in 1:nset ) {
  set = sets[[j]]
  t = table( all.genes %in% set , all.genes %in% degs )
  if( any( dim(t) != 2 ) ) next
  test = fisher.test( t , alternative="greater" )
  p[j,i] = test$p.value
  o[j,i] = t[2,2]
  or[j,i] = test$estimate
 }
 cat( "done" , i , "\n" )

 fisher.res[[k]] = list( p = p , o = o , or = or )
}
} # end main loop

names( fisher.res ) = c("bidirectional","down","up")
save( fisher.res , file="go_enrichments.fisher.limma_q05.RData")


 
p.dn = as.data.frame(fisher.res[[2]]$p)

q = apply( p.dn , 2 , p.adjust , method = "BH" )
q = as.data.frame( q )

bonf.dn = p.dn * nset * ncond
#> .05/(ncond*nset)
#[1] 6.007305e-07

setnames = names(sets)
setnames[ bonf$s129.Q111.early < 0.01 ]

colSums( bonf < 0.01 )

nsig = rowSums( q < 0.01 )
setnames[ nsig > 3 ]





