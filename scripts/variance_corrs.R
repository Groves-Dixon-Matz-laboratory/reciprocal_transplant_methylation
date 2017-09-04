#variance_corrs.R
#This script tested idea that environmentally driven increase in GBM
#correlate with increased transcriptional stability
#Found no correlation.
#Groves Dixon
#last updated 7-26-17

#### SETUP ####
library(adegenet)
setwd("~/gitreps/reciprocal_transplant_methylation/")
source("~/gitreps/reciprocal_transplant_methylation/scripts/reciprocal_methylation_project_functions.R")



##### UPLOAD DATA ####
ll=load("datasets/ge_3mo_v2.RData")
colnames(vsd) = colnames(counts)
head(vsd)
rna.vsd=vsd
samples = colnames(rna.vsd)
tks = rna.vsd[,samples[append(grep("KK", samples), grep("OK", samples))]]
tos = rna.vsd[,samples[append(grep("KO", samples), grep("OO", samples))]]
head(tks)
head(tos)
tkv = apply(tks, 1, var)
tov = apply(tos, 1, var)
vr = data.frame(tov/tkv)


ll=load("datasets/splitModels_MBD.RData")
ll=load("datasets/mbd_3mo.RData")
head(k2o.r)
head(o2k.r)
#assign plastic genes
CUT=0.01
sig.genes1 = k2o.r[k2o.r$pvalue<CUT,]
sig.genes2 = o2k.r[o2k.r$pvalue<CUT,]

sig1 = rbind(sig.genes1, sig.genes2)
head(sig1)
sig2=merge(sig1, vr, by = 0)
dim(sig2)
head(sig2)
plot(log(sig2$tov.tkv)~sig2$log2FoldChange, xlab = "Log2 GBM O vs K", ylab = "log( var(Orpheus) / var(Keppel) )")
lm1=lm(log(sig2$tov.tkv)~sig2$log2FoldChange)
summary(lm1)
abline(lm1)
upin = sig2$log2FoldChange<0
boxplot(log(sig2$tov.tkv)~upin, names = c("K", "O"), xlab="Transplant Site with Higher Methylation", ylab="log( var(Orpheus) / var(Keppel) )")








