#mbd_density_and_overlap.R
setwd("~/gitreps/reciprocal_transplant_methylation")
source("~/gitreps/reciprocal_transplant_methylation/scripts/reciprocal_methylation_project_functions.R")


#----------- grab data from full dataset tests --------------#
### origin effects ###
lnames = load('datasets/orico_GENEBODIES_meth_p-1000_200.Rdata') #file output from DEseq_gene_bodies.R.
lnames

#rename to keep separate from tag-seq results
mg.orico = orico

#note log2 fold change is transplant O vs K
head(mg.orico)

#modify column names for merging with tag-seq data
colnames(mg.orico) = paste('met', colnames(mg.orico), sep = "_")
head(mg.orico)

### transplant effects ###
lnames = load("datasets/traco_GENEBODIES_meth_p-1000_200.Rdata")
lnames
mg.traco = traco
head(mg.traco)
colnames(mg.traco) = paste('met', colnames(mg.traco), sep = "_")
head(mg.traco)



#----------- grab data from split model tests --------------#
ll=load("datasets/splitModels_MBD.RData")
head(k2o.r)   #transplant effect among keppel corals
head(o2k.r)   #transplant effect among orpheus corals
head(okatk.r) #origin effect at keppel
head(okato.r) #origin effect at orpheus


#----------- subset for tendency genes --------------#
#for full tests
CUT=0.01
ori = na.omit(data.frame(mg.orico[,c('met_log2FoldChange', 'met_pvalue')]))
trans = na.omit(data.frame(mg.traco[,c('met_log2FoldChange', 'met_pvalue')]))
sori = ori[ori$met_pvalue<CUT,]
strans = trans[trans$met_pvalue<CUT,]
head(sori)
head(trans)
nrow(sori)
nrow(strans)


#FOR SPLIT MODEL TESTS
#transplant
sig.genes1 = rownames(k2o.r)[k2o.r$pvalue<CUT];length(sig.genes1)
sig.genes2 = rownames(o2k.r)[o2k.r$pvalue<CUT];length(sig.genes2) ##increa
sig.genes = unique(append(sig.genes1, sig.genes2))
sig.genes=sig.genes[!is.na(sig.genes)]
length(sig.genes)

#for origin
sig.genes1.o = rownames(okatk.r)[okatk.r$pvalue<CUT];length(sig.genes1)
sig.genes2.o = rownames(okato.r)[okato.r$pvalue<CUT];length(sig.genes2) ##increa
sig.genes.o = unique(append(sig.genes1.o, sig.genes2.o))
sig.genes.o=sig.genes.o[!is.na(sig.genes.o)]
length(sig.genes.o)

#----------- plot mbd densities --------------#
YLIM=c(0, 0.35)
ADJ=-.25
par(mfrow=c(2,2))
plot_subset_mbd_density(rownames(sori), 'red', YLIM= YLIM, MAIN='Origin-Specific\nOO & OK vs KO & KK')
legend('topright', c('All genes (n=25733)', 'Origin-specific (n=692)'), col=c('black', 'red'), lwd=2)
mtext('A', side = 3, line = 1, adj = ADJ, cex = 2, xpd=T)

plot_subset_mbd_density(rownames(strans), 'red', YLIM= YLIM, MAIN='Environment-Specific\nOO & KO vs OK & KK')
legend('topright', c('All genes (n=25733)', 'Site-specific (n=97)'), col=c('black', 'red'), lwd=2)
mtext('B', side = 3, line = 1, adj = ADJ, cex = 2, xpd=T)

plot_subset_mbd_density(sig.genes.o,'red', YLIM= YLIM, MAIN = 'Split Origin-Specific\nOO vs KO\nKK vs OK')
legend('topright', c('All genes (n=25733)', 'Origin-specific (n=316)'), col=c('black', 'red'), lwd=2)
mtext('C', side = 3, line = 1, adj = ADJ, cex = 2, xpd=T)


plot_subset_mbd_density(sig.genes,'red', YLIM= YLIM, MAIN = 'Split Environment-Specific\nOO vs OK\nKK vs KO')
legend('topright', c('All genes (n=25733)', 'Site-specific (n=560)'), col=c('black', 'red'), lwd=2)
mtext('D', side = 3, line = 1, adj = ADJ, cex = 2, xpd=T)


#----------- overlaps --------------#
sum(rownames(sori) %in% rownames(strans))
sum(sig.genes %in% sig.genes.o)

library(VennDiagram)

plot.new()
draw.pairwise.venn(length(sig.genes), length(sig.genes.o), cross.area=sum(sig.genes %in% sig.genes.o))







