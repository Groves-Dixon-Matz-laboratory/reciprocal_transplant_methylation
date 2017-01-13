#---------------- Upload the data ------------------
#SET UP THE DATA TO RUN DESEQ
library('DESeq2')
#SAVE YOUR DIRECTORY NAME AND SET WD
directory<-"/Users/grovesdixon/lab_files/projects/recip_meth/deseq_12_13_16"
setwd(directory)

#READ IN THE COUNTS DATA. THESE ARE EXPORTED AT THE END OF THE MBD-seq data processing pipeline (see MBD-seq_Data_Processing_Walkthrough).
lnames=load("gbm_counts_p-1000_200.Rdata")
lnames
counts=gcounts

counts=counts[,order(colnames(counts))]
head(counts)
dim(counts)

#remove the flowthrough samples
mets = grep('m', colnames(counts))
counts = counts[, mets]
head(counts)
dim(counts)

#optionally remove transplants to look at origin effects
# homes = colnames(counts)[append(grep('KK', colnames(counts)), grep('OO', colnames(counts)))]
# counts = counts[, colnames(counts) %in% homes]

#LOOK AT THE MEANS OF THE COLUMNS
#HOW MANY ARE GREATER THAN 3?
dim(counts)
mns = apply(counts, 1, mean)
table(mns>10)
counts = counts[mns > 10,]
dim(counts)


#------- BUILD A DATAFRAME ASSOCIATING SAMPLE NAMES WITH TREATMENT CONDITIONS ------------
#set up sample names by clearing away .counts
sample = colnames(counts)
length(sample)


#load the trait data output from assemble_trait_data.R
traitData = read.table("assembledTraitData.tsv")
length(sample)
traitData = traitData[rownames(traitData) %in% sample,]
dim(traitData)
#double-check that names are in proper order between the datasets
sum(sample==rownames(traitData)) == length(sample)


#build the dataframe with all the varialbes
coldata <- traitData
rownames(coldata) = sample


#-------- GET VARIANCE STABILIZED COUNTS ---------
ddsHTSeq<-DESeqDataSetFromMatrix(counts,
	colData = coldata,
	design = formula(~1))
#estimate size factors
dds <- estimateSizeFactors(ddsHTSeq, type = "iterate")
dds<-estimateDispersions(dds)
dds<-nbinomWaldTest(dds)
rld = rlog(dds)


#### check the dispersion-mean rank ####
library("vsn")notAllZero <- (rowSums(counts(dds))>0)meanSdPlot(log2(counts(dds,normalized=TRUE)[notAllZero,] + 1))meanSdPlot(assay(rld[notAllZero,]))####################################



#convert to dataframe
rld.df = assay(rld)
colnames(rld.df) = colnames(counts)
#get DEseq results
res = results(dds)
#subset genes with baseMean expression above 3
res = res[res$baseMean > 3,]
# #subset the variance stabilized counts the same way
rld.df = rld.df[rownames(rld.df) %in% rownames(res),]
# save.image("iterate_local_11-30-16.Rdata")
save.image("wgcnaInitialize_p1000_200_iterate_mnCount6.Rdata")

#---------------------- IDENTIFY OUTLIERS AS IN WGCNA TUTORIAL -----------------------
directory<-"/Users/grovesdixon/lab_files/projects/recip_meth/deseq_12_13_16"
setwd(directory)
lnames=load("wgcnaInitialize_p1000_200_iterate_mnCount6.Rdata")

#heatmap
pheatmap(cor(rld.df), cluster_cols = T, cluster_rows = T, clustering_distance_rows = "maximum", clustering_distance_cols = "maximum")

#pca
#plot spidered MBD-seq pca TRANSPLANT
CEX=1.2
par(mar=c(5, 4, 4, 2) + 0.1)
ori=as.character(meth.coldata$origin)
ori[ori=="K"]<-15
ori[ori=="O"]<-19
ori=as.numeric(ori)
trans = as.character(meth.coldata$transplant)
trans[trans=="O"]<-ocol
trans[trans=="K"]<-kcol
m.pca = mod.plotPCA.df(meth.df, meth.coldata, intgroup="transplant", returnData=T,ntop= NTOP)
m.pca$PC1=m.pca$PC1*-1
m.pca$PC2=m.pca$PC2*-1
#build the plot
plot(m.pca$PC2~m.pca$PC1, col=trans, pch = 19, cex=CEX, xlab = "PC1: 32%", ylab="PC2: 6%", main = "MBDseq", axes=F)
axis(1)
axis(2, las=2)
box()
o=m.pca[,c(1,2)]
ordispider(ord=o, groups=meth.coldata$colony.id, label = F, lwd=0.5, col='black')
# ordiellipse(ord=o, groups = meth.coldata$origin, col=colors, label=F, cex=1.2, lwd=2)
# legend(x=-35, y=-30, legend=c("K>K", "K>O", "O>K", "O>O"), pch = c(15, 15, 19, 19), cex=CEX, col=c(kcol,ocol,kcol,ocol))
# legend(x=-3,y=-7, legend=c("Keppel", "Orpheus"), col=colors, pch=19, title="Transplant", cex=.75)


#=====================================================================================
#
#  Code chunk 2
# transpose the dataset you have samples as rows and genes as columns
#=====================================================================================

datExpr0 = as.data.frame(t(rld.df));
dim(datExpr0)

#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================

#check that the dataset doesn't have genes or samples with too many missing values
#these would likely represent lowly expressed genes and under sequenced samples
library(WGCNA)
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
head(gsg)


#=====================================================================================
#
#  Code chunk 4

#=====================================================================================
#removing genes that were flagged with too many missing values
#note how many genes we have right now
before = ncol(datExpr0)
print(before)


if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
     printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
     printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

nrow(datExpr0)
after = ncol(datExpr0)
print(paste(before - after, "Genes With Too Many Missing Values Were Removed"))

#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================
#now cluster samples based on gene expression to identify outliers
sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,5,2,0))
plot(sampleTree, main = "Sample Clustering", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)


#build sample heatmaps to confirm the outlier results
library(pheatmap)
quartz()
head(rld.df)
ubs = colnames(rld.df)[grep("ub", colnames(rld.df))]
rld.df = rld.df
pheatmap(cor(rld.df))
pheatmap(cor(rld.df, method = 'spearman'))
pheatmap(cor(rld.df), cluster_cols = T, cluster_rows = T, clustering_distance_rows = "maximum", clustering_distance_cols = "maximum")

save(rld.df, file='all_rld.Rdata')





#=====================================================================================
#
#  Code chunk 6
# 
#=====================================================================================

#Remove outliers by setting a branch cut threshold
# Plot a line to show the cut
cut.height = 150
abline(h = cut.height, col = "red", lty = 2);
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = cut.height, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
keepSampleNames = rownames(datExpr0)[keepSamples]
outlierNames = rownames(datExpr0)[clust==0]
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr) #number of samples left after outlier removal
print(paste(length(outlierNames), "samples were flagged as outliers and removed:"))
outlierNames
print(paste(nSamples, "samples were kept"))

#=====================================================================================

#------------ SAVE RESULTS FOR WGCNA ------------
datTraits = coldata[rownames(coldata) %in% rownames(datExpr),]
#write out normalized counts for full dataset
dim(rld.df)
out = rld.df[,!colnames(rld.df) %in% outlierNames]
dim(out)   #final input for WGCNA will have 13161 genes and 50 samples
# write.csv(out, "recipMeth_rlogWGCNA_input_baseMean3_outliersRemoved.csv")


#look at the objects to save
dim(datExpr)
dim(datTraits) 

#output data for all samples
save(datExpr, datTraits, file = "wgcna-01_output_p1000_200_iterate_mnCount6.RData") #use the Rdata file for input to run WGCNA on full dataset


#subset down to only the timepoint 2 samples
dim(datExpr)
dim(datTraits)
datExprSub=datExpr[grep("2m", rownames(datExpr)),]
datTraitsSub=datTraits[grep("2m", rownames(datTraits)),]
dim(datExprSub)
dim(datTraitsSub)
datExpr=datExprSub
datTraits=datTraitsSub
save(datExpr, datTraits, file = "wgcna-01_TimePoint2only_output.RData") #use the Rdata file for input to run WGCNA on full dataset





