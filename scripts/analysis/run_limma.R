#preliminaries
#set working directory
setwd("/proj/price1/sament/hdlux")
options(stringsAsFactors=F)
#load libraries needed in this analysis
library(limma)

#expression data
load("hdlux.full.ComBat_corrected_expr.RData")

#column metadata
load("metadata.full.RData")

#probe annotation
anno = read.delim("/proj/price1/sament/hdlux/HDLux_agilent_gene_list-2.txt")
matchProbes = match( rownames(ComBatExpr) , anno$ProbeID )
anno = anno[ matchProbes , ]

#which samples go into the current analysis?
useSamples = which(traitData2$allele_nominal %in% c("Q111" , "WT") )

expr = ComBatExpr
meta = traitData
attach( meta )
table( allele_nominal , week , sex )

n = nrow(meta)
age.3s = quantile( meta$age , c(0.33,0.67) )
age.bins = rep( NA , n )
age.bins[ meta$age <= age.3s[1] ] = "early"
age.bins[ meta$age > age.3s[1] & meta$age <= age.3s[2] ] = "middle"
age.bins[ meta$age > age.3s[2] ] = "late"


strain = meta$strain
strain = gsub( "129" , "s129" , strain )

allele = meta$allele_nominal

sex = meta$sex

group = paste( strain , allele , sex , age.bins , sep="." )

#now run a linear model using the variables that we have created
#the design matrix contains the levels of each sample for each variable
design = model.matrix( ~ 0 + group )
colnames(design) = gsub( "group" , "" , colnames(design) )

contrasts = makeContrasts(
	B6.Q111.early = (B6.Q111.F.early+B6.Q111.M.early)/2 - 
			(B6.WT.F.early+B6.WT.M.early)/2 ,
	B6.Q111.middle = (B6.Q111.F.middle+B6.Q111.M.middle)/2 - 
                        (B6.WT.F.middle+B6.WT.M.middle)/2 ,
        B6.Q111.late = (B6.Q111.F.late+B6.Q111.M.late)/2 -
                        (B6.WT.F.late+B6.WT.M.late)/2 ,
        B6.Q92.early = (B6.Q92.F.early+B6.Q92.M.early)/2 - 
                        (B6.WT.F.early+B6.WT.M.early)/2 ,
        B6.Q92.middle = (B6.Q92.F.middle+B6.Q92.M.middle)/2 -
                        (B6.WT.F.middle+B6.WT.M.middle)/2 ,
        B6.Q92.late = (B6.Q92.F.late+B6.Q92.M.late)/2 -
                        (B6.WT.F.late+B6.WT.M.late)/2 ,
        B6.Q50.early = (B6.Q50.F.early+B6.Q50.M.early)/2 - 
                        (B6.WT.F.early+B6.WT.M.early)/2 ,
        B6.Q50.middle = (B6.Q50.F.middle+B6.Q50.M.middle)/2 -
                        (B6.WT.F.middle+B6.WT.M.middle)/2 ,
        B6.Q50.late = (B6.Q50.F.late+B6.Q50.M.late)/2 -
                        (B6.WT.F.late+B6.WT.M.late)/2 ,
        CD1.Q111.early = (CD1.Q111.F.early+CD1.Q111.M.early)/2 - 
                        (CD1.WT.F.early+CD1.WT.M.early)/2 ,
        CD1.Q111.middle = (CD1.Q111.F.middle+CD1.Q111.M.middle)/2 -
                        (CD1.WT.F.middle+CD1.WT.M.middle)/2 ,
        CD1.Q111.late = (CD1.Q111.F.late+CD1.Q111.M.late)/2 -
                        (CD1.WT.F.late+CD1.WT.M.late)/2 ,
        FVB.Q111.early = (FVB.Q111.F.early+FVB.Q111.M.early)/2 - 
                        (FVB.WT.F.early+FVB.WT.M.early)/2 ,
        FVB.Q111.middle = (FVB.Q111.F.middle+FVB.Q111.M.middle)/2 -
                        (FVB.WT.F.middle+FVB.WT.M.middle)/2 ,
        FVB.Q111.late = (FVB.Q111.F.late+FVB.Q111.M.late)/2 -
                        (FVB.WT.F.late+FVB.WT.M.late)/2 ,
        s129.Q111.early = (s129.Q111.F.early+s129.Q111.M.early)/2 -
                        (s129.WT.F.early+s129.WT.M.early)/2 ,
        s129.Q111.middle = (s129.Q111.F.middle+s129.Q111.M.middle)/2 -
                        (s129.WT.F.middle+s129.WT.M.middle)/2 ,
        s129.Q111.late = (s129.Q111.F.late+s129.Q111.M.late)/2 -
                        (s129.WT.F.late+s129.WT.M.late)/2 ,
	B6.Q111 = (B6.Q111.F.early+B6.Q111.M.early+
		B6.Q111.F.middle+B6.Q111.M.middle+
		B6.Q111.F.late+B6.Q111.M.late)/6 -
		(B6.WT.F.early+B6.WT.M.early+
		B6.WT.F.middle+B6.WT.M.middle+
		B6.WT.F.late+B6.WT.M.late)/6 ,
	B6.Q92 = (B6.Q92.F.early+B6.Q92.M.early+
		B6.Q92.F.middle+B6.Q92.M.middle+
		B6.Q92.F.late+B6.Q92.M.late)/6 -
                (B6.WT.F.early+B6.WT.M.early+
                B6.WT.F.middle+B6.WT.M.middle+
                B6.WT.F.late+B6.WT.M.late)/6 ,
	B6.Q50 = (B6.Q50.F.early+B6.Q50.M.early+
		B6.Q50.F.middle+B6.Q50.M.middle+
		B6.Q50.F.late+B6.Q50.M.late)/6 -
                (B6.WT.F.early+B6.WT.M.early+
                B6.WT.F.middle+B6.WT.M.middle+
                B6.WT.F.late+B6.WT.M.late)/6 ,
	CD1.Q111 = (CD1.Q111.F.early+CD1.Q111.M.early+
		CD1.Q111.F.middle+CD1.Q111.M.middle+
		CD1.Q111.F.late+CD1.Q111.M.late)/6 -
		(CD1.WT.F.early+CD1.WT.M.early+
                CD1.WT.F.middle+CD1.WT.M.middle+
                CD1.WT.F.late+CD1.WT.M.late)/6 ,
        FVB.Q111 = (FVB.Q111.F.early+FVB.Q111.M.early+
                FVB.Q111.F.middle+FVB.Q111.M.middle+
                FVB.Q111.F.late+FVB.Q111.M.late)/6 -
                (FVB.WT.F.early+FVB.WT.M.early+
                FVB.WT.F.middle+FVB.WT.M.middle+
                FVB.WT.F.late+FVB.WT.M.late)/6 ,
        s129.Q111 = (s129.Q111.F.early+s129.Q111.M.early+
                s129.Q111.F.middle+s129.Q111.M.middle+
                s129.Q111.F.late+s129.Q111.M.late)/6 -
                (s129.WT.F.early+s129.WT.M.early+
                s129.WT.F.middle+s129.WT.M.middle+
                s129.WT.F.late+s129.WT.M.late)/6 ,
	Q111 = (B6.Q111.F.early+B6.Q111.M.early+
                B6.Q111.F.middle+B6.Q111.M.middle+
                B6.Q111.F.late+B6.Q111.M.late+
		CD1.Q111.F.early+CD1.Q111.M.early+
                CD1.Q111.F.middle+CD1.Q111.M.middle+
                CD1.Q111.F.late+CD1.Q111.M.late+
		FVB.Q111.F.early+FVB.Q111.M.early+
                FVB.Q111.F.middle+FVB.Q111.M.middle+
                FVB.Q111.F.late+FVB.Q111.M.late+
		s129.Q111.F.early+s129.Q111.M.early+
                s129.Q111.F.middle+s129.Q111.M.middle+
                s129.Q111.F.late+s129.Q111.M.late)/24 -
		(B6.WT.F.early+B6.WT.M.early+
                B6.WT.F.middle+B6.WT.M.middle+
                B6.WT.F.late+B6.WT.M.late+
		CD1.WT.F.early+CD1.WT.M.early+
                CD1.WT.F.middle+CD1.WT.M.middle+
                CD1.WT.F.late+CD1.WT.M.late+
		FVB.WT.F.early+FVB.WT.M.early+
                FVB.WT.F.middle+FVB.WT.M.middle+
                FVB.WT.F.late+FVB.WT.M.late+
		s129.WT.F.early+s129.WT.M.early+
                s129.WT.F.middle+s129.WT.M.middle+
                s129.WT.F.late+s129.WT.M.late)/24 ,
	levels=design )


# fit a linear model using limma
fit = lmFit( expr , design=design )
fit = contrasts.fit( fit , contrasts )
fit = eBayes(fit)


# extract the p-values, fold-changes, etc.
res = topTable( fit , coef = 1 , number = Inf , sort.by = "none" )
colnames(res) = paste( colnames(contrasts)[1] , colnames(res) , sep = "_" )
for( i in 2:ncol(contrasts) ) {
	top = topTable( fit , coef=i , number = Inf , sort.by = "none" )
	colnames(top) = paste( colnames(contrasts)[i] , colnames(top) , sep="_" )
	res = cbind( res , top )
	cat( i , "\n" )
}
res = cbind( anno[ , c(1,3,6) ] , res )

# a few summary statistics
fdr = grep( "adj.P.Val" , names(res) )
fdr = res[ , fdr ]
colSums( fdr < 0.05 )

pval = grep( "P.Value" , names(res) )
pval = res[ , pval ]
colSums( pval < 0.01 )

# save the results

write.table( res , file="hdlux_all_strains_limma_fit_pvals.txt" , row.names=F  , quote=F , sep="\t")
save(res , file="hdlux_all_strains_limma_fit_pvals.RData")

# intersection between contrasts
n = ncol(fdr)
o = matrix( 0 , n , n )
p = matrix( 1 , n , n )
or = matrix( 0 , n , n )
genes = res[,2]
unique_genes = unique( genes )
for( i in 1:n ) {
 for( j in 1:n ) {
  a = unique( genes[ fdr[,i] < 0.05 ] )
  b = unique( genes[ fdr[,j] < 0.05 ] )
  t = table( unique_genes %in% a , unique_genes %in% b )
  if( any(dim(t)!=2) ) next
  test = fisher.test( t , alternative = "greater" )
  o[i,j] = t[2,2]
  p[i,j] = test$p.value
  or[i,j] = test$estimate
 }
}
colnames(o) = colnames(p) = colnames(or) = rownames(o) = rownames(p) = rownames(or) = colnames(fdr)


genes[ fdr$B6.Q111_adj.P.Val < 0.05 &  fdr$B6.Q50_adj.P.Val < 0.05 ]
genes[ fdr$CD1.Q111_adj.P.Val < 0.05 &  fdr$B6.Q50_adj.P.Val < 0.05 ]
unique( genes[ fdr$FVB.Q111_adj.P.Val < 0.05 &  fdr$B6.Q50_adj.P.Val < 0.05 ])
sort( unique( genes[ fdr$B6.Q111_adj.P.Val < 0.05 &  fdr$CD1.Q111_adj.P.Val < 0.05 ]))
  

o.p01 = matrix( 0 , n , n )
p.p01 = matrix( 1 , n , n )
or.p01 = matrix( 0 , n , n )

for( in 1:n ) {
 for( j in 1:n ) {
  a = unique( genes[ fdr[,i] < 0.05 ] )
  b = unique( genes[ fdr[,j] < 0.05 ] )
  t = table( unique_genes %in% a , unique_genes %in% b )
  if( any(dim(t)!=2) ) next
  test = fisher.test( t , alternative = "greater" )
  o.p01[i,j] = t[2,2]
  p.p01[i,j] = test$p.value
  or.p01[i,j] = test$estimate
 }
}
colnames(o.p01) = colnames(p.p01) = colnames(or.p01) = rownames(o.p01) = rownames(p.p01) = rownames(or.p01) = colnames(pvals)







