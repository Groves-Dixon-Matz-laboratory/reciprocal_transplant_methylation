#DEseq_unbound_v_flowthrough.R
#Groves Dixon
#last updated 5-16-17

#Purpose: take raw MBD-seq counts, estimate absolute methylation by 
#comparing captured and flowthrough samples, save subsetted counts for 
#further analysis; plot pieces of figure 1



#---------------- Upload the data ------------------
#SET UP THE DATA TO RUN DESEQ
library('DESeq2')
#SAVE YOUR DIRECTORY NAME AND SET WD
directory<-"~/gitreps/reciprocal_transplant_methylation/"
setwd(directory)

#READ IN THE COUNTS DATA. THESE ARE EXPORTED AT THE END OF THE MBD-seq data processing pipeline (see MBD-# seq_Data_Processing_Walkthrough).
counts=read.table('datasets/promoter_gbm_counts_p-1000_200.tsv',header=TRUE,row.names=1,sep="\t"); head(counts)  #counts on just CDS regions with promoter regions excluded
counts=read.table('datasets/GENE_counts_3-6-18.tsv',header=TRUE,row.names=1,sep="\t"); head(counts)  #counts to full annotated gene regions
head(counts)
dim(counts)


#LOOK AT HOW MANY READS MAPPED TO ANNOTATED CODING SEQUENCES
wdat = counts[grep("__", rownames(counts)),]
gdat = counts[!rownames(counts) %in% rownames(wdat),]
dim(counts)
dim(wdat)
dim(gdat)

#get stats for CDS alignments only
gcounts = apply(gdat, 2, sum)
summary(gcounts)


#now remove the special cases counts
remove = rownames(counts)[grep("__", rownames(counts))]
counts = counts[!rownames(counts) %in% remove,]
dim(counts)

#gather the feature types and locus names
featureType = c()
genBank = c()
for (i in rownames(counts)){
	x=strsplit(i, "-")[[1]]
	featureType = append(featureType, x[1])
	genBank = append(genBank, x[2])
}
head(featureType)
head(genBank)
head(counts)
length(featureType)
dim(counts)

#reformat the counts data for the two genetic contexts and save

#set up data for the promoter signal
pcounts = counts[featureType == 'promoter',]
geneNames = genBank[featureType == 'promoter']
rownames(pcounts) = geneNames
head(pcounts)
dim(pcounts)
save(pcounts, file="datasets/promoter_counts_p-1000_200.Rdata")



#set up dataset for the gbm signal
gcounts = counts[featureType == 'exon',]
geneNames = genBank[featureType == 'exon']
rownames(gcounts) = geneNames
head(gcounts)
dim(gcounts)
save(gcounts, file="datasets/gbm_counts_p-1000_200.Rdata")


#now select the gene context you want to get methylation for
# counts = pcounts  #for promoters
counts = gcounts    #for gene bodies


#subset for the job1 samples that were paired ub and met samples
sample = colnames(counts)
ubs = sample[grep('ub', sample)]
mets = c()
for (i in ubs){
	mets = append(mets, paste(strsplit(i, '_', fixed = T)[[1]][1], '_2m', sep = ""))
}
keep = c(ubs,mets)
counts=counts[,colnames(counts) %in% keep]
head(counts)
dim(counts)



#LOOK AT THE MEANS OF THE COLUMNS
#HOW MANY ARE GREATER THAN 3?
mns = apply(counts, 1, mean)
table(mns > 3)

#------- BUILD A DATAFRAME ASSOCIATING SAMPLE NAMES WITH TREATMENT CONDITIONS ------------
#set up boolean variable for catpured vs flowthrough
sample=colnames(counts)
captured = sample
captured[grep('m', captured)] <- TRUE
captured[grep('ub', captured)] <- FALSE
sample
captured

colony.id = substr(sample, start=1, stop=3)


#build the dataframe with all the varialbes
coldata <- data.frame(sample, captured, colony.id)
coldata

#build DESeq input table
ddsHTSeq<-DESeqDataSetFromMatrix(counts,
	colData = coldata,
	design = formula(~colony.id + captured))

#run DESeq
dds<-estimateSizeFactors(ddsHTSeq) #decided to go iterate based on Michael Love's post that this fits a size factor that maximizes the likelihood of the null model, which is the optimal size factor "assuming there is no differential expression". Since there is little differential methylation in the dataset this seems appropriate. Since it is more flexible it will likely respond better to systemic shifts such as the entire weak-methylation component increasing.
dds <- estimateDispersions(dds, fitType="local") #local provides just slightly better correlation with CpGo/e, but it's pretty much the same
dds <- nbinomLRT(dds, reduced = ~colony.id)


resultsNames(dds)
res = results(dds, contrast=c('captured', 'TRUE', 'FALSE'))
head(res)
summary(res)
hist(res$log2FoldChange, breaks = 50)

#save the analysis
save.image('datasets/capturedVflowthroughResults_p-1000_200.Rdata')


#SAVE YOUR DIRECTORY NAME AND SET WD
# save(res, file = 'capturedVflowthroughResults.Rdata')
directory<-"~/gitreps/reciprocal_transplant_methylation/"
setwd(directory)
source("scripts/reciprocal_methylation_project_functions.R")
lnames = load('datasets/capturedVflowthroughResults_p-1000_200.Rdata')


#--------- plot figure 1 ---------------

#UPLOAD CPG DATA AND MERGE WITH MBD-SCORES
library(scales)
cdat = read.table("datasets/digitifera_rna_cpg_data.tsv", header = T) #estimated based on coding regions (see Adigitifera_annotations.txt)
cdat$cpgOE = (cdat$CpG/(cdat$C*cdat$G))*(cdat$length^2/(cdat$length-1))##equation in Gavery and Roberts (oysters)
colnames(cdat)[1]<-'genbank'
head(cdat)
#upload the rna annotation data output from adigitifera_anotation_table.py (these must be somewhere as a table but I couldn't find them)
adat = read.table("datasets/adigitifera_rna_annotations.tsv", sep = "\t", header = T) #from here: ftp://ftp.ncbi.nih.gov/genomes/Acropora_digitifera/RNA/rna.gbk.gz
adat = adat[!duplicated(adat$locusName),]
head(adat)
dim(adat)

mergedat = merge(cdat, adat, by = 'genbank')
rownames(mergedat) = mergedat$locusName


#merge with mbd-seq data based on geneID
res$genBank = rownames(res)
res2 = merge(data.frame(res), mergedat, by = 0)
head(res2)
dim(res2)

#BUILD THE PLOTS

quartz()
#plotting variables
cex.lab = 1
xlim = c(0, 1.5)
cex.axis = 1
adj = .95
alpha = 0.05

#plot distribution
par(mfrow = c(1,2))
x = hist(res$log2FoldChange, breaks = 70, xlim = c(-7, 7),main = "", xlab = 'MBD-score', mgp=MGP)
# mtext('A', side = 3, line = 0, adj = 0, cex = 1.5)

#plot correlation with CpGoe
PCH = 19
CEX = 0.2
LAS = 0
alpha = .05
plot(res2[,'cpgOE'] ~ res2[,'log2FoldChange'], xlab ="MBD-score", ylab =expression("CpG"["o/e"]), col = alpha('black', alpha), xlim=c(-7,7), ylim = c(0, 1.5), axes = F, cex.lab = cex.lab, pch = PCH, cex = CEX, las = LAS, mgp=MGP);axis(1,mgp=MGP);axis(2,mgp=MGP)
z = cor.test(res2[,'log2FoldChange'], res2[,'cpgOE'], method = "spearman")
r = signif(z$estimate, digits = 2)
title(bquote(rho ~ "=" ~ .(r[1]) * "****"), line = -1, adj = adj, font.main = 3)



#plot correlation with bisulfite data
lnames=load('datasets/methPct_mbd_score.Rdata')
#do it Mclearth 2016 style (Statistical Rethinking)
library('rethinking')
gene.mns$logm=log(gene.mns$methPct, 2)
#set up model using map
m1 <- map(
		alist(
		logm ~ dnorm(mu, sigma),
		mu <- a + b*mbd.score,
		a~dnorm(0, 10),
		b~dnorm(.1, 2),
		sigma~dunif(0,50)
		), 
		data= gene.mns)
	

precis(m1)
mbd.seq <- seq(from=-4, to=0.5, length.out=1000)
pred.dat <- list(mbd.score=mbd.seq)
mu<-link(m1, data=pred.dat)
mu.mean = apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI, prob=.95)
sim.mbd<-sim(m1, data=pred.dat)
mbd.PI<-apply(sim.mbd, 2, PI, prob=.95)

#plot results
plot(log(gene.mns$methPct,2)~gene.mns$mbd.score, ylab = expression(paste('log'[2], ' Mean % Methylation')), xlab = 'MBD-score', cex=1, axes=T,mgp=MGP)
lines(mbd.seq, mu.mean, col='red')
shade(mu.PI, mbd.seq)
shade(mbd.PI, mbd.seq)
lm2=lm(log(gene.mns$methPct,2)~gene.mns$mbd.score)
r2=sprintf("%.2f", round(summary(lm2)$r.squared, digits=2))
p=summary(lm2)$coefficients[2,4]
p
title(bquote("R"^2 ~ "=" ~ .(r2) * "***"), line = -.75, adj = .1, font.main = 3)


#----------------------------------------
#PLOT CORRELATION WITH TRANSCRIPTION




#----------- plot supplementary component figure using mclust and mixtools -----------
library(mclust)
library(Hmisc)
library(mixtools)
library(mixtools)

#grab the mbd-scores
dim(res)
res=na.omit(res)
dim(res)
x=res$log2FoldChange
#use Bayesian Information Criterion to select the optimal mixture model/number of components
bic = mclustBIC(x, G = c(1:5), modelNames = c("V"))
summary(bic)
#fit the model with two components
mod2 = Mclust(x, G = 2, modelNames = c('V'))

par(mfrow = c(1,2))
par(mar = c(5, 4, 4, 2) + 0.1)
plot(bic, G = c(1:5))
mtext('A', side = 3, line = 1, adj = -.25, cex = 2)
##NULLs
mil2mix = emmix(x,NULL,NULL,NULL,2)
summary(mil2mix)
mtext('B', side = 3, line = 1, adj = -.25, cex = 2)

#now use mclust to get classes for each gene
mod2 = Mclust(x, G = 2, modelNames = c('V'))
classes = data.frame(rownames(res), x, mod2$class)
colnames(classes) = c('gene', 'mbd.score', 'mbd.class')
head(classes)
hiWeak = max(classes$mbd.score[classes$mbd.class == 1])
lowStrong = classes[classes$class == 2,]
lowStrong = min(classes$mbd.score[classes$mbd.score > -2])
separator=mean(c(hiWeak, lowStrong))
abline(v=separator, lwd = 2)
save(classes, file="datasets/mbd_classes_p-1000_200.Rdata")
#--------------------------------------------------------------------------------------

#use Bayesian Information Criterion to select the optimal mixture model/number of components
bic = mclustBIC(x, G = c(1:5), modelNames = c("V"))
summary(bic)
#fit the model with two components
mod2 = Mclust(x, G = 2, modelNames = c('V'))
#fit the optimal model
mod2 = Mclust(x, G = 2, modelNames = c('V'))
classes = data.frame(rownames(mod2$class




# make a plot showing the distribution of MBD-scores
# This takes a long time, so keep commented out if you don't need
#and the locations/densities of the mixture components
par(mfrow = c(2, 1))
par(mar = c(0, 4, 0, 2) + .1)
hist(mdat$mbd.score, breaks = 70, xlim = c(-5.3, 10), main = "", xlab = "\n", axes = F)
at = c(-5, 0, 5, 10)
axis(2, at = c(0, 1e3, 2e3), las = 1)
par(mar = c(5, 4, 0, 2) + .1)
plot(mod2, what = "classification", col = c('red', 'green', 'purple'), axes = F, main = "\n", xlab = "\n", xlim = c(-5.3, 10))
# abline(v = at)
axis(1, at = at)
title(ylab = 'Component')
title(xlab = "MBD-score")






######## LOOK WHERE ZOOX MODULES FALL OUT IN DISTRIBUTION ########
par(mar = c(4,4,1.5, 0) + 0.1)
par(mfrow = c(1,1))
x = hist(res$log2FoldChange, breaks = 70, xlim = c(-7, 7),main = "", xlab = 'Log2 Fold Difference')
mtext('A', side = 3, line = 0, adj = 0, cex = 1.5)

res3 = na.omit(res)
plot(density(res3$log2FoldChange))


res=res3
head(na.omit(res))
blue=read.table("blue_genes.txt", header = T)
brown=read.table('brown_genes.txt', header = T)
turq=read.table("turquoise_genes.txt", header = T)

modList = list(blue, brown, turq)
sub = res[res$locusName %in% blue$x,]
lines(density(sub$log2FoldChange), col = 'blue')
sub = res[res$locusName %in% brown$x,]
lines(density(sub$log2FoldChange), col = 'brown')
sub = res[res$locusName %in% turq$x,]
lines(density(sub$log2FoldChange), col = 'turquoise')


alpha=0.25
x = hist(res$log2FoldChange, breaks = 70, xlim = c(-7, 7),main = "", xlab = 'Log2 Fold Difference')

x = hist(classes$mbd.score, breaks = 70, xlim = c(-7, 7),main = "", xlab = 'Log2 Fold Difference')
sub = res[res$locusName %in% blue$x,]
hist(sub$log2FoldChange, col = alpha('blue', alpha), add = T, breaks = 50)
sub = res[res$locusName %in% brown$x,]
hist(sub$log2FoldChange, col = alpha('brown', alpha), add = T, breaks = 50)
sub = res[res$locusName %in% turq$x,]
hist(sub$log2FoldChange, col=alpha('turquoise', alpha), add = T, breaks = 50)

head(turq$x)
subset=turq$x
plot_subset_mbd=function(subset, color, BREAKS=50){
	x = hist(classes$mbd.score, breaks = 70, xlim = c(-7, 7),main = "", xlab = 'Log2 Fold Difference')
	sub=classes[classes$gene %in% subset,]
	hist(sub$mbd.score, col=color, add=T, breaks=50)
}
plot_subset_mbd(subset, 'blue', 50)


# #check correlation between mbd-score and read counts
# readCounts = read.table("filteredReadCounts.tsv", header = F)
# colnames(readCounts) = c('sample', 'count')
# readCounts = readCounts[grep("2m", readCounts$sample),]
# readCounts = readCounts[order(readCounts$count),]
# readCounts = readCounts[readCounts$sample %in% colnames(counts),]
# head(readCounts)
# dim(readCounts)
# lowCount = readCounts$sample[1:3]
# lowCount = append(as.character(lowCount), paste(substr(as.character(lowCount), start=1, stop=3), "2ub", sep = "_"))
# highCount = readCounts$sample[10:12]
# highCount = append(as.character(highCount), paste(substr(as.character(highCount), start=1, stop=3), "2ub", sep = "_"))
# lowCount
# highCount
# all.res = res
# head(counts)
# counts = counts[,colnames(counts) %in% lowCount]
# head(counts)
# #now go back and repeat





