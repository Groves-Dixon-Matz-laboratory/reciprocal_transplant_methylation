#correlation_with_expressionV2.R
#This script compares DESeq results between the MBD-seq and tag-seq data
setwd("~/gitreps/reciprocal_transplant_methylation")
source("~/gitreps/reciprocal_transplant_methylation/scripts/reciprocal_methylation_project_functions.R")


#----------- UPLOAD MBD-SEQ RESULTS -----------#

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


#----------- UPLOAD TAG-SEQ RESULTS -----------#

### origin effects ###
lnames = load("datasets/orico_RNAseq.Rdata") #output from DEseq_Tagseq.R. Has origin effects for OO+OK vs KK+KO.
lnames

#rename the object to keep separate from MBDseq data
rnaorico = orico

#note log2 fold change is origin O vs K
head(rnaorico)

#modify coloumn names for merging
colnames(rnaorico) = paste('rna', colnames(rnaorico), sep = "_")
head(rnaorico)

### transplant effects ###
lnames = load('datasets/traco_RNAseq.Rdata')
lnames
rnatraco = traco
head(rnatraco) #note log2 fold change is transplant O vs K
colnames(rnatraco) = paste('rna', colnames(rnatraco), sep = "_")
head(rnatraco)

#----------- MERGE MBD-SEQ AND TAG-SEQ DATA -----------#

### origin data ###
x=data.frame(mg.orico)
y=data.frame(rnaorico)
odat = merge(x,y, by = 0)
head(odat)
dim(odat)

### transplant data ###
x=data.frame(mg.traco)
y=data.frame(rnatraco)
tdat = merge(x, y, by = 0)
head(tdat)
dim(tdat)


#----------- ORIGIN VOLCANO PLOT -----------#
# volcano plot
volcano_plot(data.frame(mg.traco), pcol='met_pvalue', log2col='met_log2FoldChange', fdrcol='met_padj', XLIM=c(-4,4), LEGEND='Sig. GBM FDR<0.1')
mtext('A', side = 3, line = 1, adj = -.25, cex = 2, xpd=T)

#----------- PLOT CORRELATION BETWEEN ORIGIN DIFFERENCES IN METHYLATION AND TRANSCRIPTION -----------#

#choose plotting variables
library(scales)
CUT=0.01
P.TYPE='pvalue'
alpha = 0.2
XLIM=c(-5,5)
YLIM=c(-4,4)
text.x = 2.75
text.y=-3
CEX=.75
greys = grey.colors(2)
par(mfrow=c(1,1))

#build plot for full dataset comparisons (all KK+KO samples vs all OO+OK samples)
plot(odat$rna_log2FoldChange ~ odat$met_log2FoldChange, col = alpha('black', alpha), xlim = XLIM, ylim = YLIM, xlab = expression(paste("Log"[2], ' Difference Methylation')), ylab = expression(paste("Log"[2], " Difference Transcription")), axes = F, mgp=c(2.1,1,0), cex=CEX)
axis(1); axis(2, las=2);box()
abline(h=0, v=0, lty=2, col='grey')
lm.all = lm(odat$rna_log2FoldChange ~ odat$met_log2FoldChange)
summary(lm.all)

#subset for genes that show significant variation
sig.meth = odat[odat[,paste('met', P.TYPE, sep="_")] < CUT,]
dim(sig.meth)

#overlay points for subset
points(sig.meth$rna_log2FoldChange ~ sig.meth$met_log2FoldChange, col = 'black', pch = 21, xlim = XLIM, ylim = YLIM, bg='red', cex=CEX)

#stats
lm1 = lm(sig.meth$rna_log2FoldChange ~ sig.meth$met_log2FoldChange)
abline(lm1, col = 'red')
summary(lm1)
N=nrow(na.omit(sig.meth))
R2=paste(sprintf("%.3f", round(summary(lm1)$r.squared, 3)), "****",sep='')
z = cor.test(sig.meth$rna_log2FoldChange, sig.meth$met_log2FoldChange, method="spearman")
print(z)
rho = z$estimate
p = z$p.value
text(x=text.x, y=text.y, bquote('N =' ~ .(N)), col='red')
text(x=text.x, y=text.y-.6, bquote("R"^2 ~ '=' ~ .(R2)), col='red')
legend('topleft', 'Sig. GBM (p<0.01)', pt.bg='red', col='black', pch=21, pt.cex=1.2, box.lwd=0, inset=c(0, -.2), xpd=T)
mtext('B', side = 3, line = 1, adj = -.25, cex = 2, xpd=T)



###subset for genes that show significant variation in transcription
# sig.rna = odat[odat[,paste('rna', P.TYPE, sep="_")] < CUT,]
# dim(sig.rna)

# #overlay points for subset
# plot(odat$rna_log2FoldChange ~ odat$met_log2FoldChange, col = alpha('black', alpha), xlim = XLIM, ylim = YLIM, xlab = expression(paste("Log"[2], ' Methylation')), ylab = expression(paste("Log"[2], " Transcription")), axes = F, mgp=c(2.1,1,0))
# axis(1); axis(2, las=2);box()
# abline(h=0, v=0, lty=2, col='grey')
# lm.all = lm(odat$rna_log2FoldChange ~ odat$met_log2FoldChange)
# summary(lm.all)
# points(sig.rna$rna_log2FoldChange ~ sig.rna$met_log2FoldChange, col = 'black', pch = 21, xlim = XLIM, ylim = YLIM, bg='dodgerblue')

# #stats
# lmrna = lm(sig.rna$rna_log2FoldChange ~ sig.rna$met_log2FoldChange)
# abline(lmrna, col = 'blue')
# summary(lmrna)
# N=nrow(na.omit(sig.rna))
# R2=paste(sprintf("%.3f", round(summary(lmrna)$r.squared, 3)), "****",sep='')
# z = cor.test(sig.meth$rna_log2FoldChange, sig.meth$met_log2FoldChange, method="spearman")
# print(z)
# rho = z$estimate
# p = z$p.value
# text(x=text.x, y=text.y, bquote(atop('N =' ~ .(N), "R"^2 ~ '=' ~ .(R2))), col='dodgerblue')

#now replot and overlay double-subset to include only significant variation in tag-seq data
plot(odat$rna_log2FoldChange ~ odat$met_log2FoldChange, col = alpha('black', alpha), xlim = XLIM, ylim = YLIM, xlab = expression(paste("Log"[2], ' Difference Methylation')), ylab = expression(paste("Log"[2], ' Difference Transcription')), axes = F, mgp=c(2.1,1,0), cex=CEX)
axis(1); axis(2, las=2);box()
abline(h=0, v=0, lty=2, col='grey')
sig.both = sig.meth[sig.meth[,paste('rna', P.TYPE, sep="_")] < CUT,]
points(sig.both$rna_log2FoldChange ~ sig.both$met_log2FoldChange, col = 'black', pch = 21, bg='purple', cex=CEX*1.2)
abline(lm(sig.both$rna_log2FoldChange ~ sig.both$met_log2FoldChange), col = 'purple')
lm2=lm(sig.both$rna_log2FoldChange ~ sig.both$met_log2FoldChange)
summary(lm2)
N=nrow(na.omit(sig.both))
R2=paste(sprintf("%.3f", round(summary(lm2)$r.squared, 3)), "****",sep='')
rsymbol=expression("R"^2)
z = cor.test(sig.both$rna_log2FoldChange, sig.both$met_log2FoldChange, method = "spearman")
print(z)
rro = z$estimate
p = z$p.value
text(x=text.x, y=text.y, bquote('N =' ~ .(N)), col='purple')
text(x=text.x, y=text.y-.6, bquote("R"^2 ~ '=' ~ .(R2)), col='purple')
legend('topleft', 'Sig. GBM & GE (p<0.01)', pt.bg='purple', col='black', pch=21, pt.cex=1.2, box.lwd=0, , inset=c(0, -.2), xpd=T)
mtext('C', side = 3, line = 1, adj = -.25, cex = 2, xpd=T)
summary(lm1)
# summary(lmrna)
summary(lm2)

#----------- CHECK IF COUNTS OF SIGNIFICANT GBM GENES CORRELATES WITH SIG GE COUNTS -----------#
#use fisher's exact test to test if being significant for GBM increases probablity of significance for GE
so = sig_overlap(odat, 'met_pvalue', 'rna_pvalue', 0.01, 'Row.names')


#----------- UPLOAD SPLIT MODEL RESULTS -----------#
#purpose here is to show that differences identified in one set of clones
#can predict differential transcription in the transplanted clone-mates

### MBDseq origin results ###
lnames = load("datasets/splitModels_MBD.RData")
lnames
#note log2 fold change is transplant O vs K
head(home.r)
mg.orico=home.r

#modify column names for merging with tag-seq data
colnames(mg.orico) = paste('met', colnames(mg.orico), sep = "_")
head(mg.orico)

### tagseq origin results ###
lnames=load("datasets/splitModels_GE.RData")
lnames
colnames(okatk.r) = paste('rna', colnames(okatk.r), sep = "_")
colnames(okato.r) = paste('rna', colnames(okato.r), sep = "_")
head(okatk.r) #origin effect for corals kept at Keppel
head(okato.r) #origin effect for corals kept at Orpheus


lnames=load('datasets/transr_GE.Rdata')
colnames(trans.r) = paste('rna', colnames(trans.r), sep = "_")
head(trans.r)

#----------- MERGE SPLIT MBD-SEQ AND TAG-SEQ DATA -----------#

### for corals at Keppel ###
x=data.frame(mg.orico)
y=data.frame(okatk.r)
odat.k = merge(x,y, by = 0)
head(odat.k)
dim(odat.k)

### for corals at Orpheus ###
x=data.frame(mg.orico)
y=data.frame(okato.r)
odat.o = merge(x,y, by = 0)
head(odat.o)
dim(odat.o)

### for transplants ###
x=data.frame(mg.orico)
y=data.frame(trans.r)
odat.t = merge(x,y, by = 0)
head(odat.t)
dim(odat.t)

#----------- PLOT THE CORRELATIONS FOR SPLIT MODELS -----------#

#choose plotting variables
CUT=0.01
P.TYPE='pvalue'
alpha = 0.2
XLIM=c(-2,2)
YLIM=c(-2,2)
text.x = 1.3
text.y=-1.6
greys = grey.colors(2)
par(mfrow=c(1,2))
par(mar=c(5, 6, 4, 2) + 0.1)

### set odat to either split model to plot ###
odat=odat.k;YLAB='(OK - KK)';redP='***'; purpleP="**";letter1="A";letter2="B"
odat=odat.o; YLAB='(OO - KO)';redP="**"; purpleP="*";letter1="C";letter2="D"


odat=odat.t; YLAB='(OK - KO)';redP="**"; purpleP="&";letter1="A";letter2="B"


#build plot for full dataset comparisons (all KK+KO samples vs all OO+OK samples)
odat=odat[!is.na(odat$rna_pvalue),]
odat=odat[!is.na(odat$met_pvalue),]

plot(odat$rna_log2FoldChange ~ odat$met_log2FoldChange, col = alpha('black', alpha), xlim = XLIM, ylim = YLIM, xlab='', ylab='', axes = F, main='')
axis(1); axis(2, las=2);box()
title(xlab = expression(paste("Log"[2], ' Difference Methylation')), line=2.5)
title(xlab='(OO - KK)', line=3.5)
title(ylab=expression(paste("Log"[2], " Difference Transcription")), line=3)
title(ylab=YLAB, line=4.5)
abline(h=0, v=0, lty=2, col='grey')
lm.all = lm(odat$rna_log2FoldChange ~ odat$met_log2FoldChange)
summary(lm.all)

#subset for genes that show significant variation
sig.meth = odat[odat[,paste('met', P.TYPE, sep="_")] < CUT,]
dim(sig.meth)

#overlay points for subset
points(sig.meth$rna_log2FoldChange ~ sig.meth$met_log2FoldChange, col = 'black', pch = 21, xlim = XLIM, ylim = YLIM, bg='red')

#stats
lm1 = lm(sig.meth$rna_log2FoldChange ~ sig.meth$met_log2FoldChange)
abline(lm1, col = 'red')
summary(lm1)
# N=nrow(na.omit(sig.meth))
N=nrow(sig.meth)
R2=paste(sprintf("%.3f", round(summary(lm1)$r.squared, 3)), "",sep='')
z = cor.test(sig.meth$rna_log2FoldChange, sig.meth$met_log2FoldChange, method="spearman")
print(z)
rho = z$estimate
p = z$p.value
P=sprintf("%.4f", round(summary(lm1)$coefficients[2,4], digits=4))
mtext(letter1, side = 3, line = 1, adj = -.25, cex = 2, xpd=T)
text(x=text.x, y=text.y, bquote(atop('n =' ~ .(N), "R"^2 ~ '=' ~ .(R2) * .(redP) )), col='red')


#now replot and overlay double-subset to include only significant variation in tag-seq data
plot(odat$rna_log2FoldChange ~ odat$met_log2FoldChange, col = alpha('black', alpha), xlim = XLIM, ylim = YLIM, xlab='', ylab='', axes = F, main='')
axis(1); axis(2, las=2);box()
title(xlab = expression(paste("Log"[2], ' Difference Methylation')), line=2.5)
title(xlab='(OO - KK)', line=3.5)
title(ylab=expression(paste("Log"[2], " Difference Transcription")), line=3)
title(ylab=YLAB, line=4.5)
abline(h=0, v=0, lty=2, col='grey')
sig.both = sig.meth[sig.meth[,paste('rna', P.TYPE, sep="_")] < CUT,]
points(sig.both$rna_log2FoldChange ~ sig.both$met_log2FoldChange, col = 'black', pch = 21, cex = 1, bg='purple')
abline(lm(sig.both$rna_log2FoldChange ~ sig.both$met_log2FoldChange), col = 'purple')
lm2=lm(sig.both$rna_log2FoldChange ~ sig.both$met_log2FoldChange)
summary(lm2)
N=nrow(sig.both[!is.na(sig.both$Row.names), ])
R2=paste(sprintf("%.3f", round(summary(lm2)$r.squared, 3)), "",sep='')
rsymbol=expression("R"^2)
z = cor.test(sig.both$rna_log2FoldChange, sig.both$met_log2FoldChange, method = "spearman")
print(z)
rro = z$estimate
p = z$p.value
P=sprintf("%.4f", round(summary(lm2)$coefficients[2,4], digits=4))
mtext(letter2, side = 3, line = 1, adj = -.25, cex = 2, xpd=T)
text(x=text.x, y=text.y, bquote(atop('n =' ~ .(N), "R"^2 ~ '=' ~ .(R2) ^ .(purpleP))), col='purple')
summary(lm1)
summary(lm2)

#Do fisher test
so = sig_overlap(odat, 'met_pvalue', 'rna_pvalue', 0.05, 'met_pvalue')





