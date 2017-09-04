#SET UP THE DATA TO RUN DESEQ
library('DESeq2')
#SAVE YOUR DIRECTORY NAME AND SET WD
directory<-"~/gitreps/reciprocal_transplant_methylation"
setwd(directory)


#load the data
lnames=load("datasets/promoter_counts_p-1000_200.Rdata")
head(pcounts)
dim(pcounts)


######## set up two different datasets, one with all timepoints and one with just 3-month samples (timepoint 2)
#gather data
all.samples=colnames(pcounts)
unbound=grep('ub', all.samples)
csample= all.samples[-unbound]
time=c()
id=c()
for (s in csample){
	t=strsplit(s, "_")[[1]][2]
	i=strsplit(s, "_")[[1]][1]
	time=append(time, sub('m', '', t))
	id=append(id, i)
}
table(time)

#subset
allcc = pcounts[,csample]
dim(allcc) #all samples, 44 3-month samples and 6 6-month samples
cc= allcc[,time=='2']
dim(cc)    #timepoint2 only samples
############################################################


#prepare conditions data for timepoint 2
head(cc)

sample=names(cc)
origin=substr(sample, start=1, stop=1)
transplant=substr(sample, start=2, stop=2)
num=substr(id[time=='2'], start=3, stop=5)
colony.id=paste(origin,num,sep="")
conditions=data.frame(sample,colony.id,origin,transplant)
means=apply(cc,1,mean)
table(means>2)
counts=cc[means>2,]







home=which(conditions$origin==conditions$transplant) #subset of KK and OO samples
o2k=which(conditions$origin=="O")                    #subset for samples from Orpheus
k2o=which(conditions$origin=="K")                    #subset for samples from Keppel
okato=which(conditions$transplant=="O")              #subset for samples transplanted to Orpheus (OO and OK)
okatk=which(conditions$transplant=="K")              #subset for samples transplanted to Orpheus (KK and KO)


#SET UP MODELS
#differences between home-site samples (OO and KK)
homed=DESeqDataSetFromMatrix(counts[,home],
	colData = conditions[home,], 
	design = formula(~ origin))

#test for effect of transplantsite site on corals from orpheus
o2kd=DESeqDataSetFromMatrix(counts[,o2k],
	colData = conditions[o2k,], 
	design = formula(~ colony.id+transplant))

#test for effect of transplant site on corals from keppel
k2od=DESeqDataSetFromMatrix(counts[,k2o],
	colData = conditions[k2o,], 
	design = formula(~ colony.id+transplant))

#test for effect of originating from keppel on corals transplanted to orphues
okatod=DESeqDataSetFromMatrix(counts[,okato],
	colData = conditions[okato,], 
	design = formula(~ origin))

#test for effect of originating from orphues on corals transplanted to keppel
okatkd=DESeqDataSetFromMatrix(counts[,okatk],
	colData = conditions[okatk,], 
	design = formula(~ origin))



#RUN MODELS
#always make contrast O - K
#transplant effect for corals originating from Orpheus (OO and OK)
o2kd=estimateSizeFactors(o2kd)
o2kd=estimateDispersions(o2kd)
o2kd=nbinomWaldTest(o2kd)
resultsNames(o2kd)
o2k.r=results(o2kd, contrast = c('transplant', 'O', 'K'))
summary(o2k.r)

#check importance of contrast for p-values
ok2.Tr = results(o2kd, contrast = c('transplant', 'O', 'K'))
ok2.Gr = results(o2kd, contrast = c('colony.id', 'O8', 'O9'))
summary(ok2.Tr) #weak transplant effects
summary(ok2.Gr) #note strong gneotype effects



#transplant effect for corals originating from Keppel (KK and KO)
k2od=estimateSizeFactors(k2od)
k2od=estimateDispersions(k2od)
k2od=nbinomWaldTest(k2od)
resultsNames(k2od)
k2o.r=results(k2od, contrast = c('transplant', 'O', 'K'))
summary(k2o.r)

#origin effect for corals replanted at their original sites (KK and OO)
homed=estimateSizeFactors(homed)
homed=estimateDispersions(homed)
homed=nbinomWaldTest(homed)
resultsNames(homed)
home.r=results(homed, contrast = c('origin', 'O', 'K'))
summary(home.r)


#origin effect for corals placed at orpheus (KO and OO)
okatod=estimateSizeFactors(okatod)
okatod=estimateDispersions(okatod)
okatod=nbinomWaldTest(okatod)
resultsNames(okatod)
okato.r=results(okatod, contrast = c('origin', 'O', 'K'))
summary(okato.r)

#origin effect for corals placed at keppel (KK and OK)
okatkd=estimateSizeFactors(okatkd)
okatkd=estimateDispersions(okatkd)
okatkd=nbinomWaldTest(okatkd)
resultsNames(okatkd)
okatk.r=results(okatkd, contrast = c('origin', 'O', 'K'))
summary(okatk.r)



#save the results
# save(okatk.r,okato.r,k2o.r,o2k.r,home.r,file="datasets/splitModels_MBD_PROMOTER.RData")





#--------------------- correlate with expression ---------------------#
### tagseq transplant results ###
lnames=load("datasets/splitModels_GE.RData")
lnames
colnames(k2o.r) = paste('rna', colnames(k2o.r), sep = "_")
colnames(o2k.r) = paste('rna', colnames(o2k.r), sep = "_")
rna.k2o.r = k2o.r
rna.o2k.r = o2k.r
head(rna.k2o.r) #effect of transplantation site on corals from Keppel (rna)
head(rna.o2k.r) #effect of transplantation site on corals from Orpheus (rna)

### prep mbdseq results ###
directory<-"~/gitreps/reciprocal_transplant_methylation"
setwd(directory)
lnames = load("datasets/splitModels_MBD_PROMOTER.RData")
#note log2 fold change is transplant O vs K
colnames(k2o.r) = paste('met', colnames(k2o.r), sep = "_")
colnames(o2k.r) = paste('met', colnames(o2k.r), sep = "_")
gbm.k2o.r = k2o.r
gbm.o2k.r = o2k.r
head(gbm.k2o.r) #effect of transplantation site on corals from Keppel
head(gbm.o2k.r) #effect of transplantation site on corals from Orpheus

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
CUT=0.05
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
odat=tdat.k;YLAB='(KO - KK)';redP='&'; purpleP="*";letter1="A";letter2="B"
odat=tdat.o; YLAB='(OO - OK)';redP="**"; purpleP="&";letter1="C";letter2="D"


# ############ OPTIONALLY SEPARATE BY MBD-SCORE CLASS ####################
###this did not turn up anything interesting
# lnames = load('datasets/mbd_classes_p-1000_200.Rdata')
# strong=classes[classes$mbd.class==2,'gene']
# weak=classes[classes$mbd.class==1,'gene']
# odat2=odat
# rownames(odat2) = odat2$Row.names

# #pick which to plot
# #strong
# odat=odat2[rownames(odat2) %in% strong,]
# dim(odat)

# #weak
# odat=odat2[rownames(odat2) %in% weak,]
# dim(odat)
# ####################################

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

if (redP == "&"){
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
if (purpleP == "&"){
	text(x=text.x, y=text.y, bquote(atop('n=' ~ .(N), "R"^2 * '=' * .(R2) * ""^.(purpleP))), col='purple')
} else {
	text(x=text.x, y=text.y, bquote(atop('n=' ~ .(N), "R"^2 * '=' * .(R2) * ""*.(purpleP))), col='purple')
}


summary(lm1)
summary(lm2)

#Do fisher test
so = sig_overlap(odat, 'met_pvalue', 'rna_pvalue', 0.05, 'met_pvalue')
