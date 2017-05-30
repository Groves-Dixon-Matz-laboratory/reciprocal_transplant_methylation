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
lnames = load("datasets/orico_RNAseq.Rdata") #output from DEseq_Tagseq.R. Has origin effects for KK+KO vs OO+OK.
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


### transplant GBM ori RNA ###
x=data.frame(mg.traco)
y=data.frame(rnaorico)
tdat = merge(x, y, by = 0)
head(tdat)
dim(tdat)




#----------- ORIGIN VOLCANO PLOT -----------#
# volcano plot
volcano_plot(data.frame(tdat), pcol='met_pvalue', log2col='met_log2FoldChange', fdrcol='met_padj', XLIM=c(-4,4), LEGEND='Sig. GBM FDR<0.1')
mtext('A', side = 3, line = 1, adj = -.25, cex = 2, xpd=T)

volcano_plot(data.frame(tdat), pcol='rna_pvalue', log2col='rna_log2FoldChange', fdrcol='rna_padj', XLIM=c(-4,4), LEGEND='Sig. GBM FDR<0.1')
mtext('B', side = 3, line = 1, adj = -.25, cex = 2, xpd=T)

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



odat=tdat


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
d=data.frame(sig.meth$met_pvalue, sig.meth$rna_log2FoldChange)
N=nrow(na.omit(d))
R2=sprintf("%.3f", round(summary(lm1)$r.squared, 3))
z = cor.test(sig.meth$rna_log2FoldChange, sig.meth$met_log2FoldChange, method="spearman")
print(z)
rho = z$estimate
p = z$p.value
text(x=text.x, y=text.y, bquote('N =' ~ .(N)), col='red')
text(x=text.x, y=text.y-.6, bquote("R"^2 ~ '=' ~ .(R2)^'&'), col='red')
legend('topleft', 'Sig. GBM (p<0.01)', pt.bg='red', col='black', pch=21, pt.cex=1.2, box.lwd=0, inset=c(0, -.2), xpd=T)
mtext('C', side = 3, line = 1, adj = -.25, cex = 2, xpd=T)



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
R2=sprintf("%.3f", round(summary(lm2)$r.squared, 3))
rsymbol=expression("R"^2)
z = cor.test(sig.both$rna_log2FoldChange, sig.both$met_log2FoldChange, method = "spearman")
print(z)
rro = z$estimate
p = z$p.value
text(x=text.x, y=text.y, bquote('N =' ~ .(N)), col='purple')
text(x=text.x, y=text.y-.6, bquote("R"^2 ~ '=' ~ .(R2)^'n/s'), col='purple')
legend('topleft', 'Sig. GBM & GE (p<0.01)', pt.bg='purple', col='black', pch=21, pt.cex=1.2, box.lwd=0, , inset=c(0, -.2), xpd=T)
mtext('D', side = 3, line = 1, adj = -.25, cex = 2, xpd=T)
summary(lm1)
# summary(lmrna)
summary(lm2)

#----------- CHECK IF COUNTS OF SIGNIFICANT GBM GENES CORRELATES WITH SIG GE COUNTS -----------#
#use fisher's exact test to test if being significant for GBM increases probablity of significance for GE
so = sig_overlap(odat, 'met_pvalue', 'rna_pvalue', 0.01, 'Row.names')


#----------- UPLOAD SPLIT MODEL RESULTS -----------#
#do origin based differences in transcription predict 
#transplant based differences in GBM?

### MBDseq TRANSPLANT results ###
lnames = load("datasets/splitModels_MBD.RData")
lnames
#note log2 fold change is transplant O vs K
colnames(k2o.r) = paste('met', colnames(k2o.r), sep = "_")
colnames(o2k.r) = paste('met', colnames(o2k.r), sep = "_")
gbm.k2o.r = k2o.r
gbm.o2k.r = o2k.r
head(gbm.k2o.r) #effect of transplantation site on corals from Keppel
head(gbm.o2k.r) #effect of transplantation site on corals from Orpheus

### tagseq ORIGIN results ###
lnames=load("datasets/splitModels_GE.RData")
lnames
colnames(home.r) = paste('rna', colnames(okatk.r), sep = "_")
rna.home = home.r
head(rna.home) #origin effect for corals kept at Orpheus


#----------- MERGE SPLIT MBD-SEQ AND TAG-SEQ DATA -----------#

### for corals at Keppel ###
x=data.frame(gbm.k2o.r) #effect of transplantation on corals from Keppel
y=data.frame(rna.home)   #origin effects for transcription
odat.k = merge(x,y, by = 0)
head(odat.k)
dim(odat.k)       #if K origin genes increase in methylation at K

### for corals at Orpheus ###
x=data.frame(gbm.o2k.r)
y=data.frame(rna.home)
odat.o = merge(x,y, by = 0)
head(odat.o)
dim(odat.o)

#----------- PLOT THE CORRELATIONS FOR SPLIT MODELS -----------#

#choose plotting variables
CUT=0.01
P.TYPE='pvalue'
alpha = 0.2
XLIM=c(-2,2)
YLIM=c(-1.5,1.5)
text.x = 1.2
text.y=1
greys = grey.colors(2)
par(mfrow=c(2,2))

### set odat to either split model to plot ###
odat=odat.k;MAIN='KK vs OK';redP='*'; purpleP="**";letter1="A";letter2="B"
odat=odat.o; MAIN='OO vs KO';redP="***"; purpleP="*";letter1="C";letter2="D"

#build plot for full dataset comparisons (all KK+KO samples vs all OO+OK samples)
plot(odat$rna_log2FoldChange ~ odat$met_log2FoldChange, col = alpha('black', alpha), xlim = XLIM, ylim = YLIM, xlab = expression(paste("Log"[2], ' Methylation')), ylab = expression(paste("Log"[2], " Transcription")), axes = F, main=MAIN)
axis(1); axis(2, las=2);box()
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
d=data.frame(sig.meth$met_log2FoldChange, sig.meth$rna_log2FoldChange)
N=nrow(na.omit(d))
R2=paste(sprintf("%.3f", round(summary(lm1)$r.squared, 3)), "",sep='')
z = cor.test(sig.meth$rna_log2FoldChange, sig.meth$met_log2FoldChange, method="spearman")
print(z)
rho = z$estimate
p = z$p.value
P=sprintf("%.4f", round(summary(lm1)$coefficients[2,4], digits=4))
mtext(letter1, side = 3, line = 1, adj = -.25, cex = 2, xpd=T)
text(x=text.x, y=text.y, bquote(atop('n =' ~ .(N), "R"^2 ~ '=' ~ .(R2) * .(redP) )), col='red')



#now replot and overlay double-subset to include only significant variation in tag-seq data
plot(odat$rna_log2FoldChange ~ odat$met_log2FoldChange, col = alpha('black', alpha), xlim = XLIM, ylim = YLIM, xlab = expression(paste("Log"[2], ' Methylation')), ylab = expression(paste("Log"[2], ' Transcription')), axes = F, main=MAIN)
axis(1); axis(2, las=2);box()
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

text(x=text.x, y=text.y, bquote(atop('n =' ~ .(N), "R"^2 ~ '=' ~ .(R2) * .(purpleP))), col='purple')
# text(x=text.x, y=text.y-.15, bquote("R"^2 ~ '=' ~ .(R2)), col='purple')
# text(x=text.x, y=text.y-.32, bquote("p" ~ '=' ~ .(P)), col='purple')
summary(lm1)
summary(lm2)

#Do fisher test
so = sig_overlap(odat, 'met_pvalue', 'rna_pvalue', 0.05, 'met_pvalue')




#----------- GBM CLASS-WISE CORRELATION WITH TRANSCRIPTION -----------#
lnames = load('datasets/mbd_classes_p-1000_200.Rdata')
head(classes)
strong=classes[classes$mbd.class==2,]
weak=classes[classes$mbd.class==1,]

















