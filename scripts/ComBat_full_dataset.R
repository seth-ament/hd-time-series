setwd("/proj/price1/sament/hdlux")
options(stringsAsFactors=F)
library(sva)
library(WGCNA)
allowWGCNAThreads()
library(flashClust)

# load the expression data and metadata

load("/proj/price1/sament/hdlux/data//data_full.expr_proc_32.RData")
expr = x
rm(x)

load("/proj/price1/sament/hdlux/data/metadata_full.expr_proc_32.RData")
traitData = x
rm(x)

# outlier detection by sample hierarchical clustering (all samples)
datExpr0 = t(expr)
sampleTree = flashClust(dist(datExpr0), method = "average");
strain = as.numeric( as.factor(traitData$strain ) )
allele_nominal = as.numeric( as.factor(traitData$allele_nominal))
sex = as.numeric(as.factor(traitData$sex ))
week = traitData$week
date_hyb = as.numeric(as.factor(traitData$date_hyb))
date_harvest = as.numeric(as.factor(traitData$date_harvest))
slide_id = as.numeric(as.factor(traitData$slide_id))
datTraits = data.frame( strain , allele_nominal , sex , week , date_hyb , slide_id , date_harvest)
traitColors = numbers2colors(datTraits, signed = FALSE);
setwd("/proj/price1/sament/hdlux")
pdf("hdlux_full_sampletree_with_trait_heatmap_2015-2-20.pdf")
plotDendroAndColors(sampleTree, traitColors, groupLabels = names(datTraits), main = "Sample dendrogram and trait heatmap")
dev.off()

# some batch effects by date_hyb are apparent. run combat.

batch = traitData$date_hyb
mod = model.matrix(~ as.factor( traitData$strain ) + as.factor( traitData$allele_nominal ) + as.factor( traitData$week ) )
ComBatExpr = ComBat(dat = t(datExpr0) , batch = batch , mod = mod , par.prior = TRUE , prior.plots = F )

# re-plot hierarchical clustering after combat correction.

sampleTree2 = flashClust(dist(t(ComBatExpr)), method = "average");
strain = as.numeric( as.factor(traitData2$strain ) )
allele_nominal = as.numeric( as.factor(traitData2$allele_nominal))
sex = as.numeric(as.factor(traitData2$sex ))
week = traitData2$week
date_hyb = as.numeric(as.factor(traitData2$date_hyb))
date_harvest = as.numeric(as.factor(traitData2$date_harvest))
slide_id = as.numeric(as.factor(traitData2$slide_id))
datTraits2 = data.frame( allele_nominal , sex , week , date_hyb , slide_id , date_harvest)
traitColors2 = numbers2colors(datTraits2, signed = FALSE);
pdf("hdlux_combat_fuil_sampletree_with_trait_heatmap.pdf")
par(cex = 0.6);
par(mar = c(0,4,2,0))
plotDendroAndColors(sampleTree2, traitColors2, groupLabels = names(datTraits2), main = " ComBat-corrected, no-outlier-removed samples and traits")
dev.off()

# batch effects are no longer apparent.

save(ComBatExpr , file="hdlux.full.ComBat_corrected_expr.RData")
save(traitData , file="metadata.full.RData")

