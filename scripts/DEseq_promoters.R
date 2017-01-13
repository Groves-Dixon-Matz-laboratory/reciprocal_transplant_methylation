#---------------- Upload the data ------------------
#SET UP THE DATA TO RUN DESEQ
library('DESeq2')
#SAVE YOUR DIRECTORY NAME AND SET WD
directory<-"/Users/grovesdixon/lab_files/projects/recip_meth/deseq_12_13_16"
setwd(directory)

#READ IN THE COUNTS DATA. THESE ARE EXPORTED AT THE END OF THE MBD-seq data processing pipeline (see MBD-seq_Data_Processing_Walkthrough).
lnames=load("promoter_counts_p500_a0.Rdata")
lnames=load("promoter_counts_p-500.Rdata")
lnames=load("promoter_counts_p-1000_200.Rdata")
counts=pcounts
head(counts)
dim(counts)

#remove the flowthrough and timepoint 3 samples
mets = grep('2m', colnames(counts))
counts = counts[, mets]
head(counts)
dim(counts)

#optionally remove transplants to look at origin effects
# homes = colnames(counts)[append(grep('KK', colnames(counts)), grep('OO', colnames(counts)))]
# counts = counts[, colnames(counts) %in% homes]

#LOOK AT THE MEANS OF THE COLUMNS
#HOW MANY ARE GREATER THAN 3?
mns = apply(counts, 1, mean)
# counts=counts[mns>1,] #get rid of most promoters that show little or no methylation
table(mns > 3)
dim(counts)

#---------------- set up dataframe for sample varibles -----------------------
#BUILD A DATAFRAME ASSOCIATING SAMPLE NAMESWITH TREATMENT CONDITIONS
#set up sample names by clearing away .counts
sample = c()
for (i in colnames(counts)){
	sample = append(sample, strsplit(i, '.', fixed = T)[[1]][1])
}
sample

#set up boolean variable for catpured vs flowthrough
captured = sample
captured[grep('2m', captured)] <- TRUE
captured[grep('2ub', captured)] <- FALSE
sample
captured

set up tranplant variable
transplant = sample
transplant[grep('KK', transplant)] <- 'K'
transplant[grep('OK', transplant)] <- 'K'
transplant[grep('KO', transplant)] <- 'O'
transplant[grep('OO', transplant)] <- 'O'
transplant

#set up origin variable
origin = sample
origin[grep('KK', origin)] <- 'K'
origin[grep('OK', origin)] <- 'O'
origin[grep('KO', origin)] <- 'K'
origin[grep('OO', origin)] <- 'O'
origin

#put in colony.id
tag = sub("_2m", "", sample)
tag1 = substr(tag, start=1, stop=1)
tag2 = substr(tag, start=3, stop = 5)
colony.id=paste(tag1,tag2,sep="")
colony.id

#build the dataframe with all the varialbes
COLDATA <- data.frame(sample, colony.id, transplant,origin)
COLDATA

##### FIRST TEST FOR EFFECT OF TRANSPLANT #####
#here we can account for colony.id in the model
dds<-DESeqDataSetFromMatrix(counts,
	colData = COLDATA, 
	design = formula(~ colony.id + transplant))

#test for effect of transplant controlling for colony.id
dds <- estimateSizeFactors(dds)#note found that iterate works better. See wgcna1_initialize.R type="iterate"
# save(ddsco, file="iterated_sizeFactors_PROMOTERS_colonyid-transplant_ddsco_DESeq_environment.Rdata")
dds <- estimateDispersions(dds)
dds<-nbinomLRT(dds, reduced=~colony.id)
resultsNames(dds)
traco = results(dds, contrast = c('transplant', 'O', 'K'), independentFiltering=T)
traco=traco[order(traco$pvalue),]
head(traco)
summary(traco)
save(traco, file="traco_meth_PROMOTERS_transplant-colonyid_p250both.Rdata")

rld=rlog(dds)
rld.df = assay(rld)
colnames(rld.df) = colnames(counts)
save(rld.df, file="rld_PROMOTERS__p-1000_200.Rdata")

### compare read counts for select genes
#single gene
g="LOC107337099"
d=plotCounts(ddsco, gene=g, intgroup="transplant", returnData=F)
d$sample=COLDATA$sample
d

#all significant
x=na.omit(traco)
sig=x[x$padj<0.1,]
named_plotcounts(ddsco, sig, INTGROUP='transplant')

rld=rlog(dds)


######### LRT testing #########
dds<-DESeqDataSetFromMatrix(counts,
	colData = COLDATA, 
	design = formula(~ transplant+origin))

#prepare deseq object then split for two analyses
dds <- estimateSizeFactors(dds, type = 'iterate')
dds <- estimateDispersions(dds)
dds.t<-dds #for transplant
dds.o<-dds #for origin

#lrt test for transplant
dds.t<-nbinomLRT(dds.t, reduced=~origin)
resultsNames(dds.t)
traco=results(dds.t, contrast = c('transplant', 'O', 'K'))
summary(traco)   #no differentially methylated genes for this test
head(na.omit(traco))
save(traco, file="traco_PROMOTERS_meth_p-1000_200.Rdata")

#lrt test for origin
dds.o<-nbinomLRT(dds.o, reduced=~transplant)
resultsNames(dds.o)
orico=results(dds.o, contrast = c('origin', 'O', 'K'), independentFiltering=F)
orico=orico[order(orico$pvalue),]
summary(orico)
head(na.omit(orico))
save(orico, file="orico_PROMOTERS_meth_p-1000_200.Rdata")

#### plot the read Counts ###
#all significant
source("../reciprocal_methylation_project_functions.R")
x=na.omit(orico)
sig=x[x$padj<0.1,];nrow(sig)
named_plotcounts_multipanel(dds.o, sig, INTGROUP='origin', LINE=T)

#save the origin object
save(orico, file="orico_meth_p-500.Rdata")
save(traco, file="traco_meth.Rdata")


# testing simulated
sdsco.o=DESeq(ddsimco,test="LRT",reduced=~transplant)
sdsco.t=DESeq(ddsimco,test="LRT",reduced=~origin)
sorico=results(sdsco.o)
straco=results(sdsco.t)
summary(sorico)
summary(straco)
#as expected these don't show anything

#now we need to do empiricle FDR
#first for the origin results
fdroco=fdrTable(orico$pvalue,sorico$pvalue)
quartz()
par(mfrow=c(1,2))
fdrBiCurve(fdroco,main="origin")
efdr=empiricalFDR(fdroco,plot=T,main="origin")
mtext(paste("10% FDR at p =",signif(efdr,2)),cex=0.8)
table(orico$pval<efdr) # empirical FDR based DEGs 159
orico$efdr=(orico$pval<efdr)
orico$efdr[is.na(orico$efdr)]=FALSE
head(orico)




#now repeat for the transplant results
fdrtco=fdrTable(traco$pvalue,straco$pvalue)
quartz()
par(mfrow=c(1,2))
fdrBiCurve(fdrtco,main="transplant")
efdr=empiricalFDR(fdrtco,plot=T,main="transplant")
mtext(paste("10% FDR at p =",signif(efdr,2)),cex=0.8)
table(traco$pval<efdr) # empirical FDR based DEGs 159
traco$efdr=(traco$pval<efdr)
traco$efdr[is.na(traco$efdr)]=FALSE


save(orico, file="orico_meth.Rdata")
save(traco, file='taco_meth.Rdata')

#save/load
# save.image("postEmpiricalFDR.Rdata")
directory<-"/Users/grovesdixon/lab_files/projects/recip_meth/deseq_11-4-16/dig_combo_mapped"
setwd(directory)
lnames=load("postEmpiricalFDR.Rdata")

###############################
####### output for GO #########
dat = orico

sign = dat[,'log2FoldChange'] > 0
sign[sign == TRUE] <- 1
sign[sign == FALSE] <- -1
out = na.omit(data.frame(rownames(dat), dat$stat*sign))
colnames(out) = c('locusName', 'stat')
head(out)
dim(out)
pdat = read.table("prot_annotations.tsv", header = T)
dim(pdat)
pdat = pdat[unique(pdat$locusName),] #reduce multiple protein entries from same gene to single line
head(pdat)
dim(pdat)
dim(out)
out = merge(out, pdat, by = "locusName", all=F)
dim(out)
out = out[,colnames(out) %in% c('genbank', 'stat')]
out = out[,order(colnames(out))]
head(out)
dim(out)
out=na.omit(out)
write.csv(out, 'origin_GO_input.csv', row.names = F, quote = F)



###############################
# principal [blank] analysis etc

rl=rlog(ddsco)
vsdrl=assay(rl)
colnames(vsdrl)=names(countsco)
row.names(vsdrl)=row.names(countsco)

gt=sub("_2m","",COLDATA$sample)
gt=sub("KK|KO","K",gt)
gt=sub("OO|OK","O",gt)
COLDATA$gt=gt

#save/load
# save(ddsco,vsdrl,rl,COLDATA,orico,traco,file="captureOnly.RData")
source('/Users/grovesdixon/lab_files/coding4people/lafire/DESeq_all_counts_8-17-16/multivariate_functions.R')
directory<-"/Users/grovesdixon/lab_files/projects/recip_meth/deseq_7-29-16/map_to_coding_seqs"
setwd(directory)
lnames=load("captureOnly.RData")
lnames




library(DESeq2)
library(ggplot2)
NTOP = 25000
plotPCA(rld,intgroup=c("oriK","transK"),ntop=25000)+theme_bw()
plotPCA(rl,intgroup=c("origin"),ntop= 25000)+theme_bw()







#check if any PCs separate origin
for (i in 2:8){
	mod.plotPCA(rl, intgroup=c('origin'), ntop = NTOP, returnData = T, pcs = 10, pc1=1, pc2 = i, main = 'Origin')
}
for (i in 2:8){
	mod.plotPCA(rl, intgroup=c('transplant'), ntop = NTOP, returnData = T, pcs = 10, pc1=1, pc2 = i, main = 'transplant')
}

#build sample heatmaps
library(pheatmap)
colnames(vsdrl)=paste(COLDATA$gt,COLDATA$transplant,sep=">")
pheatmap(cor(vsdrl))
pheatmap(cor(vsdrl, method = 'spearman'))


library(ape)
library(vegan)

# mns=apply(vsdrl,1,mean)
# q50=quantile(mns,prob=0.5)
# sds=apply(vsdrl,1,sd)
# q50sd=quantile(sds,prob=0.5)

co.pcoa=pcoa(vegdist(t(vsdrl),method="manhattan")/1000)
scores=co.pcoa$vectors

quartz()
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(scores[,1], scores[,2],col=as.numeric(as.factor(COLDATA$origin)),pch=19+as.numeric(as.factor(COLDATA$transplant)),mgp=c(2.1,1,0),xlab="PCo1",ylab="PCo2")
ordispider(scores,COLDATA$gt,label=F)
legend("bottomright",inset=c(-0.3,0), pch=c(20,21),legend=c("tr:K","tr:O"),bty="n")
legend("right",inset=c(-0.33,0),pch=20,col=c("black","red"),legend=c("ori:K","ori:O"),bty="n")

library(rgl)
plot3d(scores[,1], scores[,2], scores[,3],col=as.numeric(as.factor(COLDATA$origin)),type="s",radius=0.25*sqrt(as.numeric(as.factor(COLDATA$transplant))),xlab="", ylab="",zlab="")

adonis(t(vsdrl)~origin*transplant,data=COLDATA,method="manhattan")
                  # Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
# origin             1  37046784 37046784 1.58691 0.14526  0.063 .
# transplant         1  20781421 20781421 0.89018 0.08148  0.582  
# origin:transplant  1  10443864 10443864 0.44737 0.04095  0.982  
# Residuals          8 186762036 23345254         0.73230         
# Total             11 255034105                  1.00000         


#-----------------------
# PCA

mns=apply(vsdrl,1,mean)
q50=quantile(mns,prob=0.5)
sds=apply(vsdrl,1,sd)
q50sd=quantile(sds,prob=0.8)

scores = mod.plotPCA.df(vsdrl, COLDATA, intgroup = 'origin')
co.pca=prcomp(t(vsdrl))
summary(co.pca)
scores=co.pca$rotation
varexp=co.pca$sdev^2/sum(co.pca$sdev^2)

quartz()
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(scores[,1], scores[,2],col=as.numeric(as.factor(COLDATA$origin)),pch=19+as.numeric(as.factor(COLDATA$transplant)),mgp=c(2.1,1,0),xlab="PC1",ylab="PC2")
ordispider(scores,COLDATA$gt,label=T)
legend("bottomright",inset=c(-0.3,0), pch=c(20,21),legend=c("tr:K","tr:O"),bty="n")
legend("right",inset=c(-0.33,0),pch=20,col=c("black","red"),legend=c("ori:K","ori:O"),bty="n")

plot(scores[,1], scores[,2],col=as.numeric(as.factor(COLDATA$origin)),pch=as.character(COLDATA$transplant),cex=0.8)
ordispider(scores,gt,label=F)
legend("bottomright",pch=20,col=c("black","red"),legend=c("ori:K","ori:O"),bty="n")

library(rgl)
plot3d(scores[,1], scores[,2], scores[,3],col=as.numeric(as.factor(COLDATA$origin)),type="s",radius=0.04*as.numeric(as.factor(COLDATA$transplant)),xlab="", ylab="",zlab="")


#try plotting other components
plot(scores[,1], scores[,9],col=as.numeric(as.factor(COLDATA$origin)), pch = 19,mgp=c(2.1,1,0),xlab="PC1",ylab="PC2")




#SET UP THE COLUMN DATA FOR THE DESeq TABLE
# colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c("O","K"))

#----------------------

# trying without O3 and K1


countsco=countsco[,!(COLDATA$gt %in% c("K1","O3"))]
COLDATA=COLDATA[!(COLDATA$gt %in% c("K1","O3")),]
COLDATA
names(countsco)

#BUILD A DESeq INPUT TABLE FROM THE DATAFRAME
ddsco<-DESeqDataSetFromMatrix(countsco,
	colData = COLDATA, 
	design = formula(~ transplant+origin))

# generating simulated counts under null model
dds0co<-DESeqDataSetFromMatrix(countsco,
	colData = COLDATA, 
	design = ~1)
dds0co=DESeq(dds0co)
library(empiricalFDR.DESeq2)
simco=simulateCounts(dds0co)
ddsimco=DESeqDataSetFromMatrix(assay(simco),
	colData = COLDATA, 
	design = formula(~ transplant+origin))

# LRT testing
ddsco.o=DESeq(ddsco,test="LRT",reduced=~transplant)
ddsco.t=DESeq(ddsco,test="LRT",reduced=~origin)
orico=results(ddsco.o)
traco=results(ddsco.t)
summary(orico)
summary(traco)

# testing simulated
sdsco.o=DESeq(ddsimco,test="LRT",reduced=~transplant)
sdsco.t=DESeq(ddsimco,test="LRT",reduced=~origin)
sorico=results(sdsco.o)
straco=results(sdsco.t)
summary(sorico)
summary(straco)

fdroco=fdrTable(orico$pvalue,sorico$pvalue)
quartz()
par(mfrow=c(1,2))
fdrBiCurve(fdroco,main="origin")
efdr=empiricalFDR(fdroco,plot=T,main="origin")
mtext(paste("10% FDR at p =",signif(efdr,2)),cex=0.8)
table(orico$pval<efdr) # empirical FDR based DEGs 159
orico$efdr=(orico$pval<efdr)
orico$efdr[is.na(orico$efdr)]=FALSE

fdrtco=fdrTable(traco$pvalue,straco$pvalue)
quartz()
par(mfrow=c(1,2))
fdrBiCurve(fdrtco,main="origin")
efdr=empiricalFDR(fdrtco,plot=T,main="origin")
mtext(paste("10% FDR at p =",signif(efdr,2)),cex=0.8)
table(traco$pval<efdr) # empirical FDR based DEGs 159
orico$efdr=(orico$pval<efdr)
orico$efdr[is.na(orico$efdr)]=FALSE

# plotting log fold changes with or without flow-through
quartz()
par(mfrow=c(1,3))
plot(orico[,"log2FoldChange"]~ori[,"log2FoldChange"],col=rgb(0,0,0,0.1),pch=19,cex=0.7,mgp=c(2.1,1,0),xlab="capture+flow-through",ylab="capture only",main="log2 fold change\n~origin")
abline(0,1,col="red")
plotMA(ddsco.t)
plotMA(ddsco.o)

# volcano plots
quartz()
plot(-log(pvalue)~log2FoldChange,orico,col="white",mgp=c(2.1,1,0),main="~origin")
points(-log(pvalue)~log2FoldChange,subset(orico, !efdr),pch=18,col=rgb(0,0,0,0.3),cex=0.7)
points(-log(pvalue)~log2FoldChange,subset(orico, efdr),pch=18,col=rgb(1,0,0,0.3),cex=0.7)
quartz()
plot(-log(pvalue)~log2FoldChange,traco,pch=18,col=rgb(0,0,0,0.3),cex=0.7,ylim=c(0,70),mgp=c(2.1,1,0),main="~transplant")







