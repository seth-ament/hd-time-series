# create a figure showing the expression of Mlh1 across strains and genotypes

load("hdlux.full.ComBat_corrected_expr.RData")

load("metadata.full.RData")
meta = traitData


n = nrow(meta)
age.3s = quantile( meta$age , c(0.33,0.67) )
age.bins = rep( NA , n )
age.bins[ meta$age <= age.3s[1] ] = "early"
age.bins[ meta$age > age.3s[1] & meta$age <= age.3s[2] ] = "middle"
age.bins[ meta$age > age.3s[2] ] = "late"

strain = meta$strain

allele = meta$allele_nominal

sex = meta$sex

group = paste( strain , allele , sep="." )

anno = read.delim("HDLux_agilent_gene_list-2.txt")
matchProbes = match( rownames(ComBatExpr) , anno$ProbeID )
anno = anno[ matchProbes , ]


e = ComBatExpr

Mlh1 = which( anno$GeneSymbol == "Mlh1" )

Mlh1 = ComBatExpr[Mlh1,]

require( doBy )

df = data.frame( strain, allele , sex , age.bins , group , Mlh1 )
exclude = which( group == "CD1.Q50" )
df = df[ -exclude , ]

g = as.factor( df$group )
g = factor( g , c("B6.WT","B6.Q50","B6.Q92","B6.Q111"," ",
  "CD1.WT","CD1.Q111","  " ,
  "FVB.WT","FVB.Q111", "   " ,
  "129.WT","129.Q111" ) )
  
pdf("figures/Mlh1.boxplot.pdf" , height = 5 , width = 7 )
par( las = 3 , bty = "l" , cex = 1 , cex.axis = 1.5 , mar = c(7,5,5,1) )
boxplot( df$Mlh1 ~ g )
mtext( side = 3 , line = 1 , "Mlh1" , font = 4 , adj = 0 , las = 0 , cex = 2  )
mtext( side = 2 , line = 3 , "log2( Expression )" , cex = 1.5 )
dev.off()





