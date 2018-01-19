
#### SETUP ####
library(adegenet)
setwd("~/gitreps/reciprocal_transplant_methylation/")
source("~/gitreps/reciprocal_transplant_methylation/scripts/reciprocal_methylation_project_functions.R")

##### UPLOAD DATA ####
#DEseq results for GE
ll=load("datasets/splitModels_GE.RData")
ll=load("datasets/ge_3mo.RData")
home.ge = data.frame(home.r)
k2o.ge = data.frame(k2o.r)
o2k.ge = data.frame(o2k.r)

colnames(home.ge) = paste("ge", colnames(home.r), sep='.')
colnames(k2o.ge) = paste("ge", colnames(k2o.r), sep='.')
colnames(o2k.ge) = paste("ge", colnames(o2k.r), sep='.')
head(home.ge)
head(k2o.ge)
head(o2k.ge)

#DEseq results for GBM
ll=load("datasets/splitModels_MBD.RData")
ll=load("datasets/mbd_3mo.RData")
head(k2o.r)
head(o2k.r)


#assign plastic genes
CUT=0.01
sig.genes1 = rownames(k2o.r)[k2o.r$pvalue<CUT];length(sig.genes1)
sig.genes2 = rownames(o2k.r)[o2k.r$pvalue<CUT];length(sig.genes2) ##increa
sig.genes = unique(append(sig.genes1, sig.genes2))
sig.genes=sig.genes[!is.na(sig.genes)]
length(sig.genes)


#CORRELATE TRANSPLANTATION EFFECT ON GBM WITH 
#for KO gbm changes
s0 = data.frame(k2o.r)[rownames(k2o.r) %in% sig.genes1,];head(s0);dim(s0)
s1 = home.ge[rownames(home.ge) %in% sig.genes1,];head(s1);dim(s1)
m=merge(s0, s1, by = 0);head(m);dim(m)
head(m)
plot(m$log2FoldChange~m$ge.log2FoldChange)

#for OK gbm changes
s0 = data.frame(o2k.r)[sig.genes2,];head(s0);dim(s0)
s1 = home.ge[sig.genes2,];head(s1);dim(s1)
m=merge(s0, s1, by = 0);head(m);dim(m)
plot(m$log2FoldChange~m$ge.log2FoldChange)


#CORRELATE TRANSPLANTATION EFFECT ON GBM WITH EFFECT ON GE
#for k2o
head(k2o)












#CORRELATE TRANSPLANTATION EFFECT ON GBM WITH 
#for KO gbm changes
s0 = data.frame(k2o.r)[rownames(k2o.r) %in% sig.genes1,];head(s0);dim(s0)
s1 = k2o.ge[rownames(k2o.ge) %in% sig.genes1,];head(s1);dim(s1)
m=merge(s0, s1, by = 0);head(m);dim(m)
plot(m$log2FoldChange~m$ge.log2FoldChange)

#for OK gbm changes
s0 = data.frame(o2k.r[sig.genes2,])
s1 = home.ge[sig.genes2,]
m=merge(s0, s1, by = 0)
plot(m$log2FoldChange~m$ge.log2FoldChange)









m=merge(k2o.r, k2o.ge, by = 0)
plot(m$log2FoldChange~m$ge.log2FoldChange)
lm1=lm(m$log2FoldChange~m$ge.log2FoldChange)
abline(lm1, col='red')
