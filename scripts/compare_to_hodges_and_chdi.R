# compare DEGs in time series study to Hodges and CHDI Allelic Series

hom = read.delim("/proj/price1/sament/resources/HOM_MouseHumanSequence.rpt")
human = hom[ hom$Common.Organism.Name == "human" , ]
mouse = hom[ hom$Common.Organism.Name == "mouse, laboratory" , ]
orth = merge( human[,c(1,4)] , mouse[,c(1,4)] , by=1 )
colnames(orth)[2:3] = c("human","mouse")
orth = orth[,2:3]

# Hodges (human post mortem caudate)
load("/proj/price1/CHDI/users/sament/allelic_series_striatum/independent_datasets/human/hodges/trn_2015-8-2/limma_pvals.RData")
hodges = res
hodges = hodges[ order( hodges$P.Value ) , ]
hodges = hodges[ duplicated(hodges[,2]) == F , ]
hodges = hodges[,c(2,5,8,9)]
hodges = merge( orth , hodges , by = 1 )

# CHDI allelic series (knock-in mice)
load("/proj/price1/CHDI/users/sament/allelic_series_striatum/edgeR.lrt.RData")
chdi = res
chdi = chdi[,c(1,39,42,43)]

load("../hdlux_all_strains_limma_fit_pvals.RData")
hdlux = res
hdlux = hdlux[,c(2,16,19,20)]
hdlux = hdlux[ order( hdlux[,3] ) , ]
hdlux = hdlux[ duplicated(hdlux[,1]) == F , ]

# HDLux B6.Q111 (late) vs. CHDI B6.Q111 6mo
y = merge( chdi , hdlux , by = 1 )
t = table( y$FDR.q111.6m < 0.05 , y$B6.Q111.late_adj.P.Val < 0.05 )
fisher.test( t , alternative="greater" )$p.value

hdlux = merge( orth , hdlux , by.x = 2 , by.y = 1 )
chdi = merge( orth , chdi , by.x = 2 , by.y = 1 )

x = merge( chdi , hdlux , by=1:2 )
x = merge( x , hodges , by.x=2:1 , by.y=1:2 )

hodges.q05 = which( x$adj.P.Val < 0.05 )
hdlux.q05 = which( x$B6.Q111.late_adj.P.Val < 0.05 )
chdi.q05 = which( x$FDR.q111.6m < 0.05 )

degs = matrix( 0 , nrow=nrow(x) , ncol=3 )
degs[ hdlux.q05 , 1 ] = 1
degs[ chdi.q05 , 2 ] = 1
degs[ hodges.q05 , 3 ] = 1
colnames(degs) = c("b6.q111.late","chdi.q111.6m","hodges")

library(limma)
vennCounts(degs)
pdf("venn.hodges.chdi.pdf")
vennDiagram( degs )
dev.off()
 
# HDLux vs. Hodges
t = table( degs[,1] == 1 , degs[,3] == 1 )
fisher.test( t , alternative="greater" )$p.value
# p = 3.1e-14


load("../hdlux_all_strains_limma_fit_pvals.RData")
hdlux = res

genes = setdiff( hdlux[,2] , c("",NA))

b6.q111.old = setdiff(unique( hdlux[ hdlux$B6.Q111.late_adj.P.Val < 0.05 ,2]) ,c(NA,"") )

b6.q111.mid =setdiff(unique(hdlux[ hdlux$B6.Q111.middle_adj.P.Val < 0.05 ,2]) ,c(NA,"") )

cd1.q111.old =setdiff(unique( hdlux[ hdlux$CD1.Q111.late_adj.P.Val < 0.05 ,2]) ,c(NA,"") )

degs = matrix( 0 , nrow=length(genes) , ncol=3 )
degs[ genes %in% b6.q111.old , 1 ] = 1
degs[ genes %in% b6.q111.mid , 2 ] = 1
degs[ genes %in% cd1.q111.old , 3 ] = 1

t = table( degs[,1] == 1 , degs[,2] == 1 )
fisher.test( t , alternative="greater" )$p.value
# p = 1.0e-174

t = table( degs[,1] == 1 , degs[,3] == 1 )
fisher.test( t , alternative="greater" )$p.value
# p = 1.4e-173

t = table( degs[,2] == 1 , degs[,3] == 1 )
fisher.test( t , alternative="greater" )$p.value
# 2.0e-94

(121+77) / (121+77+81)
(77+101) / (77+101+44)











