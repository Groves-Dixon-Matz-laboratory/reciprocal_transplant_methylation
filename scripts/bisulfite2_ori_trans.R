#bisulfite2_ori_trans.R
#Build plots and run statistics to validate origin and transplant effects from MBD-seq data
#See bisulfite_seq_metadata.excel for information on the loci that were tested


setwd("~/gitreps/reciprocal_transplant_methylation/bisulfite_validation/")
source("~/gitreps/reciprocal_transplant_methylation/scripts/bisulfite_functions.R")
source("~/gitreps/reciprocal_transplant_methylation/scripts/reciprocal_methylation_project_functions.R")
lnames = load("../datasets/ProteinTable.Rdata")
lnames = load("../datasets/deseqObjects_GENEBODIES_promoter1000_200.Rdata")

library(DESeq2)
library(reshape)
library(plotrix)


#---------- LOOK AT GENES WITH EXPECTED ORIGIN EFFECTS ----------#
lnames=load("../datasets/bisulfite_dat.Rdata")
lnames
#locus 13
par(mfrow=c(1,2))
res=plot.mean.meth(d, 'locus13', 'ori', 'methPct', 'Origin')
gene.meth.barplot(d, 'locus13', 'ori', 'By Origin')


#locus 14
res=plot.mean.meth(d, 'locus14', 'ori', 'methPct', 'Origin')
gene.meth.barplot(d, 'locus14', 'ori', 'By Origin', plus=2)


#locus 15 (only has a single Orpheus sample, do not use)
# res=plot.mean.meth(dat, 'locus15', 'ori', 'methPct', 'Origin')
# gene.meth.barplot(dat, 'locus15', 'ori', 'By Origin')


#locus 16
res=plot.mean.meth(d, 'locus16', 'ori', 'methPct', 'Origin')
gene.meth.barplot(d, 'locus16', 'ori', 'By Origin', plus=.75)

#locus 17 (only has a single representative originating from Orpheus do not use)
# res=plot.mean.meth(d, 'locus17', 'ori', 'methPct', 'Origin')
# gene.meth.barplot(d, 'locus17', 'ori', 'By Origin')

#locus 19
res=plot.mean.meth(d, 'locus19', 'ori', 'methPct', 'Origin')
gene.meth.barplot(d, 'locus19', 'ori', 'By Origin', plus=.75)

#plot the mbd-seq data for each of the above genes
bsgenes
ori.genes = c('LOC107327073', 'LOC107356898', 'LOC107356899', 'LOC107347512')
for (g in ori.genes){
	mdat = named_plotcounts(dds.o, na.omit(orico)[g,], INTGROUP='origin')
}



#---------- LOOK AT GENES WITH EXPECTED TRANSPLANT EFFECTS ----------#

#locus18
##Note that locus 18 only worked for samples originating in Keppels.
#Don't trust this one
# res=plot.mean.meth(dat, 'locus18', 'trans', 'methPct', 'Transplant')
# plot.mean.pair.diff(dat, 'locus18', 'trans', 'methPct', 'Origin', 'K')
# gene.meth.barplot(dat, 'locus18', 'trans', 'Transplant')

#locus23
par(mfrow=c(1,3))
res=plot.mean.meth(d, 'locus23', 'trans', 'methPct', 'Transplant')
plot.mean.pair.diff(d, 'locus23', 'trans', 'methPct', 'Origin', 'K')
g=gene.meth.barplot(d, 'locus23', 'trans', 'Transplant', plus=.75)



#locus24
res=plot.mean.meth(d, 'locus24', 'trans', 'methPct', 'Transplant')
plot.mean.pair.diff(d, 'locus24', 'trans', 'methPct', 'Origin', 'K')
gene.meth.barplot(d, 'locus24', 'trans', 'Transplant')

#locus25
res=plot.mean.meth(d, 'locus25', 'trans', 'methPct', 'Transplant')
plot.mean.pair.diff(d, 'locus25', 'trans', 'methPct', 'Origin', 'K')
gene.meth.barplot(d, 'locus25', 'trans', 'Transplant', plus=.1)

#plot the mbd-seq data for these
bsgenes
trans.genes = bsgenes[bsgenes$effect == 'transplant', 'gene']
for (g in trans.genes){
	mdat = named_plotcounts(dds.t, na.omit(traco)[g,], INTGROUP='transplant')
}

#make a legend
plot.new()
legend("center", c("Orpheus", "Keppel"), fill=c('dodgerblue', 'firebrick'), title="Origin")
plot.new()
legend("center", c("Orpheus", "Keppel"), fill=c('dodgerblue', 'firebrick'), title="Transplantation Site")

#---------- TEST CORRELATION BETWEEN DIFFERENCE IN MEAN MBD-SEQ COUNT AND BISULFITE METH PCT ----------#

#### FOR ORIGIN EFFECTS ###

#set up two different gene ids
all.genes = bsgenes$gene
bs.names = bsgenes$primerset
ori.names = bsgenes$primerset[bsgenes$effect=='origin']

#get differences for MBDseq and bisulfite seq for each gene
mbd.diffs = c()
bs.diffs = c()
for (g in all.genes){
	#get difference of means for mbd normalized counts
	mdat = named_plotcounts(dds.o, na.omit(orico)[g,], INTGROUP='origin')
	mns = tapply(mdat$count, mdat$origin, mean)
	diff = mns['K'] - mns['O']
	mbd.diffs = append(mbd.diffs, diff)
	#get the difference in mean percent methylation
	bs.locus=bsgenes[bsgenes$gene==g,'locus']
	lsub = d[d$locus == bs.locus,]
	mns = tapply(lsub$methPct, lsub$ori, mean)
	diff = mns['K'] - mns['O']
	bs.diffs = append(bs.diffs, diff)
}

#assemble results
ores = data.frame(mbd.diffs, bs.diffs, all.genes, bs.names)
colors = ores$bs.names
colors[colors %in% ori.names]<-19
colors[!colors %in% ori.names]<-1
colnames(ores) = c('mbd.diffs', 'bs.diffs', 'gene', 'bs.names')
lm1=lm(ores$mbd.diffs~ores$bs.diffs)
R2=round(summary(lm1)$r.squared, digits=2)
summary(lm1)

#plot
plot(ores$mbd.diffs~ores$bs.diffs, xlab = "Mean(K) - Mean(O) Bisulfite", ylab = "Mean(K) - Mean(O) MBD-seq", pch=colors, cex=2)
abline(lm1, col='red')
legend('bottomright', c('significant', 'not significant'), title='Origin Effect', pt.bg=c('black', 'white'), pch=21)
text(x=-4.8, y=28, bquote("R"^2 ~ "=" ~ .(R2) * "***"), cex=1.2)

#### FOR TRANSPLANT EFFECTS ###
all.genes = bsgenes$gene
bs.names = bsgenes$primerset
trans.names = bsgenes$primerset[bsgenes$effect=='transplant']


#get differences for MBDseq and bisulfite seq for each gene
mbd.diffs = c()
bs.diffs = c()
for (g in all.genes){
	#get difference of means for mbd normalized counts
	mdat = named_plotcounts(dds.t, na.omit(traco)[g,], INTGROUP='transplant')
	mns = tapply(mdat$count, mdat[,'transplant'], mean)
	diff = mns['K'] - mns['O']
	mbd.diffs = append(mbd.diffs, diff)
	#get the difference in mean percent methylation
	bs.locus=bsgenes[bsgenes$gene==g,'locus']
	lsub = d[d$locus == bs.locus,]
	mns = tapply(lsub$methPct, lsub$trans, mean)
	diff = mns['K'] - mns['O']
	bs.diffs = append(bs.diffs, diff)
}

#assemble results
tres = data.frame(mbd.diffs, bs.diffs, all.genes, bs.names)
colors = tres$bs.names
colors[!colors %in% trans.names]<-1
colors[colors %in% trans.names]<-19
colnames(tres) = c('mbd.diffs', 'bs.diffs', 'gene', 'bs.names')
lm1=lm(tres$mbd.diffs~tres$bs.diffs)
R2=round(summary(lm1)$r.squared, digits=2)
summary(lm1)

#plot
plot(tres$mbd.diffs~tres$bs.diffs, xlab = "Mean(K) - Mean(O) Bisulfite", ylab = "Mean(K) - Mean(O) MBD-seq", pch=colors, cex=2)
abline(lm1, col='red')
legend('topleft', c('significant', 'not significant'), title='Transplant Effect', pt.bg=c('black', 'white'), pch=21)
text(x=-3, y=7, bquote("R"^2 ~ "=" ~ .(R2)^"n/s"), cex=1.2)




#assemble combined
cres = rbind(na.omit(ores), na.omit(tres))
hist(cres$mbd.diffs)
hist(cres$bs.diffs)
lm1=lm(cres$mbd.diffs~cres$bs.diffs)
summary(lm1)
R2=round(summary(lm1)$r.squared, digits=2)


#plot
plot(cres$mbd.diffs~cres$bs.diffs, xlab = "Mean(K) - Mean(O) Bisulfite", ylab = "Mean(K) - Mean(O) MBD-seq", cex=2)
abline(lm1, col='red')
text(x=-5, y=25, bquote("R"^2 ~ "=" ~ .(R2) * "****"), cex=1.2)

#jacknife
head(cres)
par(mfrow=c(3,5))
plot(cres$mbd.diffs~cres$bs.diffs, xlab = "Mean(K) - Mean(O) Bisulfite", ylab = "Mean(K) - Mean(O) MBD-seq", cex=1.5)
lm1=lm(cres$mbd.diffs~cres$bs.diffs)
p=round(summary(lm1)$coefficients[2,4], digits=3)
text(x=-2,y=-60, labels=paste('p=',p))
summary(lm1)
abline(lm1, col='red')
title(main='All Bisulfite Loci')
title(main=paste("\n\nR2=", round(summary(lm1)$r.squared, digits=3)))
slopes=c(summary(lm1)$coefficients[2,1])
for (i in unique(cres$bs.names)){
	sub=cres[cres$bs.names != i,]
	plot(sub$mbd.diffs~ sub$bs.diffs, xlab = "Mean(K) - Mean(O) Bisulfite", ylab = "Mean(K) - Mean(O) MBD-seq", cex=1.5)
	lm1=lm(sub$mbd.diffs~ sub$bs.diffs)
	print(summary(lm1))
	abline(lm1, col='red')
	slp=summary(lm1)$coefficients[2,1]
	p=round(summary(lm1)$coefficients[2,4], digits=5)
	title(main=paste(i, 'missing'))
	title(main=paste("\n\nR2=", round(summary(lm1)$r.squared, digits=3)))
	# text(x=-2,y=-60, labels=paste('p=',p))
	legend('bottomright', legend=bquote('p = '*.(p)))
	slopes=append(slopes, slp)
}
slopes
par(mfrow=c(1,1))
hist(slopes)
