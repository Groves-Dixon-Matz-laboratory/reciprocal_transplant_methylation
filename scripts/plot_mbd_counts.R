#plot_mbd_counts.R
#code for plotting normalized counts for GBM by transplant and origin


library(DESeq2)
setwd("~/git_Repositories/reciprocal_transplant_methylation/")

#--------- load functions and data ------------------------------
source("scripts/reciprocal_methylation_project_functions.R")
lnames = load("datasets/ProteinTable.Rdata")
lnames = load("datasets/deseqObjects_GENEBODIES_promoter1000_200.Rdata")
lnames
#----------------------------------------------------------------


#----------- plot significant gene sets -------------------------
#look at significant genes
x=na.omit(traco) #for transplant
x=na.omit(orico) #for origin

#subset
sig=x[x$pvalue<0.004,];nrow(sig)

#plot significant 
named_plotcounts_multipanel(dds.o, sig, INTGROUP='origin', LINE=T)
named_plotcounts_multipanel(dds.t, sig, INTGROUP='transplant', LINE=T)
#----------------------------------------------------------------

#----------- plot individual gene -------------------------------
g="LOC107334334"
named_plotcounts(dds.t, na.omit(traco)[g,], INTGROUP='transplant')
named_plotcounts(dds.o, na.omit(orico)[g,], INTGROUP='origin')
#----------------------------------------------------------------

#----------- plot for bisulfite tested genes --------------------
gdat=read.csv("datasets/bisulfite_genes.csv");head(gdat)
#plot for transplant effect
for (i in 1:length(gdat$gene)){
	g=gdat$gene[i]
	pset=gdat$primerset[i]
	named_plotcounts(dds.t, na.omit(traco)[g,], INTGROUP='transplant')
	title(paste('\n\n\n\nprimer set', pset))
}
#plot for origin effect
for (i in 1:length(gdat$gene)){
	g=gdat$gene[i]
	pset=gdat$primerset[i]
	named_plotcounts(dds.t, na.omit(orico)[g,], INTGROUP='origin')
	title(paste('\n\n\n\nprimer set', pset))
}
#----------------------------------------------------------------




#----------- plot for any set of genes --------------------------
lnames=load("~/Desktop/turquoise.Rdata")
subset = tcors$locusName[1:30]

#for origin
x=orico[rownames(orico) %in% subset,]
named_plotcounts_multipanel(dds.o, x, INTGROUP='origin', LINE=T)

#for transplant
x=traco[rownames(traco) %in% subset,]
named_plotcounts_multipanel(dds.t, x, INTGROUP='transplant', LINE=T)

################################################