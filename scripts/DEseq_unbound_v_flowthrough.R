#---------------- Upload the data ------------------
#SET UP THE DATA TO RUN DESEQ
library('DESeq2')
#SAVE YOUR DIRECTORY NAME AND SET WD
directory<-"/Users/grovesdixon/lab_files/projects/recip_meth/deseq_12_13_16"
setwd(directory)

#READ IN THE COUNTS DATA. THESE ARE EXPORTED AT THE END OF THE MBD-seq data processing pipeline (see MBD-seq_Data_Processing_Walkthrough).
counts=read.table('promoter_gbm_counts_p500_a0.tsv',header=TRUE,row.names=1); head(counts)  #this one has promoter assigned as 1000 bp overlapping the tss
counts=read.table('promoter_gbm_counts_p-500.tsv',header=TRUE,row.names=1); head(counts)  #this one has promoter assigned as position -1000 through -500 to leave buffer from gbm signal
counts=read.table('promoter_gbm_counts_p900-100.tsv',header=TRUE,row.names=1,sep="\t"); head(counts)  #this one has promoter assigned as position -1000 through -500 to leave buffer from gbm signal
counts=read.table('promoter_gbm_counts_p250_both.tsv',header=TRUE,row.names=1,sep="\t"); head(counts)  #small promoter boundary just 250 to both sides of promoter

counts=read.table('promoter_gbm_counts_p-1000_200.tsv',header=TRUE,row.names=1,sep="\t"); head(counts)  #small promo
head(counts)
dim(counts)
remove = rownames(counts)[grep("__", rownames(counts))]
counts = counts[!rownames(counts) %in% remove,]
dim(counts)

#gather the feature types and genbank names
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
save(pcounts, file="promoter_counts_p-1000_200.Rdata")



#set up dataset for the gbm signal
gcounts = counts[featureType == 'exon',]
geneNames = genBank[featureType == 'exon']
rownames(gcounts) = geneNames
head(gcounts)
dim(gcounts)
save(gcounts, file="gbm_counts_p-1000_200.Rdata")


#now select the gene context you want to get methylation for
counts = pcounts
counts = gcounts


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
save.image('capturedVflowthroughResults_p-1000_200.Rdata')

#SAVE YOUR DIRECTORY NAME AND SET WD
directory<-"/Users/grovesdixon/lab_files/projects/recip_meth/deseq_12_13_16"
setwd(directory)
lnames = load('capturedVflowthroughResults_p-1000_200.Rdata')
lnames
# all.res = res
# low.res = res
# high.res = res
# head(low.res)
# head(all.res)
# head(high.res)
# par(mfrow=c(3,1))
# plot(density(na.omit(high.res$log2FoldChange)), col = 'red', main = 'density mbd scores')
# lines(density(na.omit(low.res$log2FoldChange)), col = 'blue')
# lines(density(na.omit(all.res$log2FoldChange)), col = 'black')
# legend(x=5,y=.2, legend=c('high', 'low', 'all'), fill=c('red', 'blue', 'black'))


#--------- plot figure 1 ---------------
#upload the iso2seq table which we'll need later
# save(res, file = 'capturedVflowthroughResults.Rdata')
directory<-"/Users/grovesdixon/lab_files/projects/recip_meth/deseq_11-4-16/dig_combo_mapped"
setwd(directory)
source("../../reciprocal_methylation_project_functions.R")
# load('capturedVflowthroughResults.Rdata')


#UPLOAD CPG DATA AND MERGE WITH MBD-SCORES
cdat = read.table("digitifera_rna_cpg_data.tsv", header = T) #estimated based on coding regions (Amil_CD)
cdat$cpgOE = (cdat$CpG/(cdat$C*cdat$G))*(cdat$length^2/(cdat$length-1))##equation in Gavery and Roberts (oysters)
colnames(cdat)[1]<-'genbank'
head(cdat)
#upload the rna annotation data output from adigitifera_anotation_table.py (these must be somewhere as a table but I couldn't find them)
adat = read.table("rna_annotations.tsv", sep = "\t", header = T)
adat = adat[!duplicated(adat$locusName),]
head(adat)
dim(adat)

mergedat = merge(cdat, adat, by = 'genbank')
rownames(mergedat) = mergedat$locusName

#upload protein annotations
# protdat = read.table("prot_annotations.tsv", sep = "\t", header = T)
# protdat = protdat[!duplicated(protdat$locusName),]
# head(protdat)
# mergedat = merge(protdat, adat, by = 'locusName')
# colnames(mergedat) = c('locusName', 'genBank', 'gi', 'rna.genBank', 'rna.gi')
# head(mergedat)

#merge with mbd-seq data based on geneID
res$genBank = rownames(res)
res2 = merge(data.frame(res), mergedat, by = 0)
head(res2)
dim(res2)
# colnames(cdat)[1]='rna.genBank'
# res3 = merge(res2, cdat, by = 'rna.genBank')
# head(res3)
# dim(res3)

#BUILD THE PLOTS

quartz()
#plotting variables
cex.lab = 1
xlim = c(0, 1.5)
cex.axis = 1
adj = .95
alpha = 0.05

#plot distribution
par(mar = c(4,4,1.5, 0) + 0.1)
par(mfrow = c(1,2))
x = hist(res$log2FoldChange, breaks = 70, xlim = c(-7, 7),main = "", xlab = 'Log2 Fold Difference')
mtext('A', side = 3, line = 0, adj = 0, cex = 1.5)

#plot correlation with CpGoe
PCH = 19
CEX = 0.2
LAS = 0
alpha = .1
par(mar = c(4,3.5,1.5,1) + 0.1)
r = plot.lm('log2FoldChange', 'cpgOE', res2, '\n', '\n', limits = T, c(-7, 7), ylim = c(0, 1.5), point.color = 'black', print.line = F)
title(paste(paste('r =', r[1]), "", sep = ""), line = -1, adj = adj, font.main = 3)
title(ylab = expression("CpG"["o/e"]), line = 2.25, cex.lab = cex.lab)
title(xlab = 'MBD-score', cex.lab = cex.lab, line = 2.5)
mtext('B', side = 3, line = 0, adj = 0, cex = 1.5)
#----------------------------------------




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
save(classes, file="/Users/grovesdixon/lab_files/projects/recip_meth/deseq_12_13_16/mbd_classes_p-1000_200.Rdata")
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
and the locations/densities of the mixture components
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





