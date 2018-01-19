#correlation_with_expression.R
#This script compares DESeq results between the MBD-seq and tag-seq data
setwd("~/gitreps/reciprocal_transplant_methylation")
source("~/gitreps/reciprocal_transplant_methylation/scripts/reciprocal_methylation_project_functions.R")


#----------- UPLOAD MBD-SEQ RESULTS -----------#

### transplant effects ###
lnames = load("datasets/traco_GENEBODIES_meth_p-1000_200.Rdata")
lnames
mg.traco = traco
head(mg.traco)
colnames(mg.traco) = paste('met', colnames(mg.traco), sep = "_")
head(mg.traco)


#----------- UPLOAD TAG-SEQ RESULTS -----------#

### transplant effects ###
lnames = load('datasets/traco_RNAseq.Rdata')
lnames
rnatraco = traco
head(rnatraco) #note log2 fold change is transplant O vs K
colnames(rnatraco) = paste('rna', colnames(rnatraco), sep = "_")
head(rnatraco)

#----------- MERGE MBD-SEQ AND TAG-SEQ DATA -----------#

### transplant data ###
x=data.frame(mg.traco)
y=data.frame(rnatraco)
tdat = merge(x, y, by = 0)
head(tdat)
dim(tdat)


#----------- TRANSPLANT VOLCANO PLOT -----------#
# volcano plot
volcano_plot(data.frame(mg.traco), pcol='met_pvalue', log2col='met_log2FoldChange', fdrcol='met_padj', XLIM=c(-4,4), LEGEND='Sig. GBM FDR<0.1')
mtext('A', side = 3, line = 1, adj = -.25, cex = 2, xpd=T)

#plot density to go with it
x=mg.traco$met_log2FoldChange
plot(density(na.omit(x)), xlim=c(-4,4), axes=F, main='',xlab='');axis(1)


#----------- PLOT CORRELATION BETWEEN ORIGIN DIFFERENCES IN METHYLATION AND TRANSCRIPTION -----------#

#choose plotting variables
library(scales)
CUT=0.01
P.TYPE='pvalue'
alpha = 0.2
XLIM=c(-5,5)
YLIM=c(-4,4)
text.x = 2.75
text.y=3.75
CEX=.75
letter.adj = -.25
greys = grey.colors(2)
par(mfrow=c(1,1))


#switch names
odat = tdat


odat = tdat[rownames(tdat) %in% rownames(strong),]
odat = tdat[rownames(tdat) %in% rownames(weak),]


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
N=nrow(na.omit(sig.meth[,c('met_log2FoldChange', 'rna_log2FoldChange')]))
R2=sprintf("%.3f", round(summary(lm1)$r.squared, 3))
z = cor.test(sig.meth$rna_log2FoldChange, sig.meth$met_log2FoldChange, method="spearman")
print(z)
rho = z$estimate
p = sprintf("%.3f", round(z$p.value, digits=3))
text(x=text.x, y=text.y, bquote('N =' ~ .(N)), col='red')
text(x=text.x, y=text.y-.6, bquote("R"^2 ~ '=' ~ .(R2)^"n/s"), col='red')
legend('topleft', 'Sig. GBM (p<0.01)', pt.bg='red', col='black', pch=21, pt.cex=1.2, box.lwd=0, inset=c(0, -.2), xpd=T)
mtext('B', side = 3, line = 1, adj = letter.adj, cex = 2, xpd=T)

#build density plots
plot(density(na.omit(odat$met_log2FoldChange)), main='', xlab='', xlim=XLIM, axes=F, col='black');axis(1)
lines(density(na.omit(sig.meth$met_log2FoldChange)), col='red')


#now replot and overlay double-subset to include only significant variation in tag-seq data
plot(odat$rna_log2FoldChange ~ odat$met_log2FoldChange, col = alpha('black', alpha), xlim = XLIM, ylim = YLIM, xlab = expression(paste("Log"[2], ' Difference Methylation')), ylab = expression(paste("Log"[2], ' Difference Transcription')), axes = F, mgp=c(2.1,1,0), cex=CEX)
axis(1); axis(2, las=2);box()
abline(h=0, v=0, lty=2, col='grey')
sig.both = sig.meth[sig.meth[,paste('rna', P.TYPE, sep="_")] < CUT,]
points(sig.both$rna_log2FoldChange ~ sig.both$met_log2FoldChange, col = 'black', pch = 21, bg='purple', cex=CEX*1.2)
abline(lm(sig.both$rna_log2FoldChange ~ sig.both$met_log2FoldChange), col = 'purple')
lm2=lm(sig.both$rna_log2FoldChange ~ sig.both$met_log2FoldChange)
summary(lm2)
N=nrow(na.omit(sig.both[,c('met_log2FoldChange', 'rna_log2FoldChange')]))
R2=sprintf("%.3f", round(summary(lm2)$r.squared, 3))
rsymbol=expression("R"^2)
z = cor.test(sig.both$rna_log2FoldChange, sig.both$met_log2FoldChange, method = "spearman")
print(z)
rro = z$estimate
p = z$p.value
text(x=text.x, y=text.y, bquote('N =' ~ .(N)), col='purple')
text(x=text.x, y=text.y-.6, bquote("R"^2 ~ '=' ~ .(R2)^"n/s"), col='purple')
legend('topleft', 'Sig. GBM & GE (p<0.01)', pt.bg='purple', col='black', pch=21, pt.cex=1.2, box.lwd=0, , inset=c(0, -.2), xpd=T)
mtext('C', side = 3, line = 1, adj = letter.adj, cex = 2, xpd=T)
summary(lm1)
# summary(lmrna)
summary(lm2)


#build density plots
plot(density(na.omit(sig.both$met_log2FoldChange)), main='', xlab='', xlim=XLIM, axes=F, col='purple');axis(1)
lines(density(na.omit(odat$met_log2FoldChange)), col='black')

#----------- CHECK IF COUNTS OF SIGNIFICANT GBM GENES CORRELATES WITH SIG GE COUNTS -----------#
#use fisher's exact test to test if being significant for GBM increases probablity of significance for GE
so = sig_overlap(odat, 'met_pvalue', 'rna_pvalue', 0.01, 'Row.names')




#----------- GBM CLASS-WISE CORRELATION WITH TRANSCRIPTION -----------#
lnames = load('datasets/mbd_classes_p-1000_200.Rdata')
head(classes)
strong=classes[classes$mbd.class==2,]
weak=classes[classes$mbd.class==1,]









#----------- UPLOAD SPLIT MODEL RESULTS -----------#
#purpose here is to show that differences identified in one set of clones
#can predict differential transcription in the transplanted clone-mates

### MBDseq transplant results ###
lnames = load("datasets/splitModels_MBD.RData")
# lnames = load("datasets/splitModels_MBD_PROMOTER.RData") #!OPTIONALLY RUN SAME PLOTS WITH PROMOTER COUNTS (build these plots in split_model_MBD_PROMOTERS.R)
lnames
#note log2 fold change is transplant O vs K
colnames(k2o.r) = paste('met', colnames(k2o.r), sep = "_")
colnames(o2k.r) = paste('met', colnames(o2k.r), sep = "_")
gbm.k2o.r = k2o.r
gbm.o2k.r = o2k.r
head(gbm.k2o.r) #effect of transplantation site on corals from Keppel
head(gbm.o2k.r) #effect of transplantation site on corals from Orpheus


### tagseq transplant results ###
lnames=load("datasets/splitModels_GE.RData")
lnames
colnames(k2o.r) = paste('rna', colnames(k2o.r), sep = "_")
colnames(o2k.r) = paste('rna', colnames(o2k.r), sep = "_")
rna.k2o.r = k2o.r
rna.o2k.r = o2k.r
head(rna.k2o.r) #effect of transplantation site on corals from Keppel (rna)
head(rna.o2k.r) #effect of transplantation site on corals from Orpheus (rna)


#----------- MERGE SPLIT MBD-SEQ AND TAG-SEQ DATA -----------#

### for corals from Keppel ###
x=data.frame(gbm.k2o.r)
y=data.frame(rna.k2o.r)
tdat.k = merge(x,y, by = 0)
head(tdat.k)
dim(tdat.k)

### for corals at Orpheus ###
x=data.frame(gbm.o2k.r)
y=data.frame(rna.o2k.r)
tdat.o = merge(x,y, by = 0)
head(tdat.o)
dim(tdat.o)

#----------- PLOT THE CORRELATIONS FOR SPLIT MODELS -----------#

#choose plotting variables
CUT=0.01
P.TYPE='pvalue'
alpha = 0.2
XLIM=c(-2,2)
YLIM=c(-1.5,1.5)
text.x = 1.1
text.y=1.2
greys = grey.colors(2)
par(mfrow=c(2,2))
par(mar=c(5, 6, 4, 2) + 0.1)

### set odat to either split model to plot ###
odat=tdat.k;YLAB='(KO - KK)';redP='n/s'; purpleP="n/s";letter1="A";letter2="B"
odat=tdat.o; YLAB='(OO - OK)';redP="****"; purpleP="n/s";letter1="C";letter2="D"

#build plot for full dataset comparisons (all KK+KO samples vs all OO+OK samples)
plot(odat$rna_log2FoldChange ~ odat$met_log2FoldChange, col = alpha('black', alpha), xlim = XLIM, ylim = YLIM, xlab='', ylab='', axes = F, main='')
axis(1); axis(2, las=2);box()
title(xlab = expression(paste("Log"[2], ' Difference Methylation')), line=2.5)
title(xlab=YLAB, line=3.5)
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
N=nrow(sig.meth[!is.na(sig.meth $Row.names), ])
R2=paste(sprintf("%.3f", round(summary(lm1)$r.squared, 3)), "",sep='')
z = cor.test(sig.meth$rna_log2FoldChange, sig.meth$met_log2FoldChange, method="spearman")
print(z)
rho = z$estimate
p = z$p.value
P=sprintf("%.4f", round(summary(lm1)$coefficients[2,4], digits=4))
mtext(letter1, side = 3, line = 1, adj = -.25, cex = 2, xpd=T)

if (redP == "n/s"){
	text(x=text.x, y=text.y, bquote(atop('n=' * .(N), "R"^2 * '=' * .(R2) * ''^.(redP) )), col='red')
} else {
	text(x=text.x, y=text.y, bquote(atop('n=' * .(N), "R"^2 * '=' * .(R2) * ''*.(redP) )), col='red')
}

#now replot and overlay double-subset to include only significant variation in tag-seq data
plot(odat$rna_log2FoldChange ~ odat$met_log2FoldChange, col = alpha('black', alpha), xlim = XLIM, ylim = YLIM, xlab='', ylab='', axes = F, main='')
axis(1); axis(2, las=2);box()
title(xlab = expression(paste("Log"[2], ' Difference Methylation')), line=2.5)
title(xlab=YLAB, line=3.5)
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

text(x=text.x, y=text.y, bquote(atop('n=' ~ .(N), "R"^2 * '=' * .(R2) * ""^.(purpleP))), col='purple')
summary(lm1)
summary(lm2)

#Do fisher test
so = sig_overlap(odat, 'met_pvalue', 'rna_pvalue', 0.05, 'met_pvalue')


#-------------------- LOOK MORE CLOSELY AT OO VS OK DIFFERENCES --------------------#

### set odat to either split model to plot ###
odat=tdat.k;MAIN='KO vs KK';redP='n/s'; purpleP="n/s";letter1="A";letter2="B"
odat=tdat.o; MAIN='OO vs OK';redP="****"; purpleP="n/s";letter1="C";letter2="D"

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



head(sig.meth)
upK = sig.meth[sig.meth$met_log2FoldChange < 0,]
upO = sig.meth[sig.meth$met_log2FoldChange > 0,]
par(mfrow=c(2,1))
hist(upK$rna_log2FoldChange, xlab='RNA log2 Fold Change');abline(v=mean(upK$rna_log2FoldChange), lwd=2, lty=2)
hist(upO$rna_log2FoldChange,xlab='RNA log2 Fold Change');abline(v=mean(upO$rna_log2FoldChange), lwd=2, lty=2)



