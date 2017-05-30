#DEseq_gene_bodies.R
#This script uses DESeq2 to perform likelihood ratio tests for
#effects of origin and transplant on GBM. These tests include
#all samples in each comparison
#transplant test:  KK + OK vs OO + KO
#origin test:      KK + KO vs OO + OK

#origin 

#---------------- Upload the data ------------------
#SET UP THE DATA TO RUN DESEQ
library('DESeq2')
#SAVE YOUR DIRECTORY NAME AND SET WD
directory<-"~/gitreps/reciprocal_transplant_methylation"
setwd(directory)

#READ IN THE COUNTS DATA. THESE ARE EXPORTED AT THE END OF THE MBD-seq data processing pipeline (see MBD-seq_Data_Processing_Walkthrough).
lnames=load("datasets/gbm_counts_p-1000_200.Rdata")

#check data
lnames
counts=gcounts
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

#set up tranplant variable
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
x = c()
for (i in sample){
	x = append(x, strsplit(i, "_")[[1]][1])
}
y=substr(x, start=1, stop=1)
z=substr(x, start=3, stop=5)
colony.id = paste(y,z,sep="")
length(unique(colony.id))

#put in treatment
treat=substr(sample, start=1, stop=2)


#build the dataframe with all the varialbes
COLDATA <- data.frame(sample, colony.id, transplant,origin,treat)
COLDATA
# save(COLDATA, counts, file='raw_counts_sample_data.Rdata')


######### LRT testing #########
dds<-DESeqDataSetFromMatrix(counts,
	colData = COLDATA, 
	design = formula(~ transplant+origin))

meth.rld=rlog(dds)
meth.coldata = COLDATA



#prepare deseq object then split for two analyses
dds <- estimateSizeFactors(dds, type = "iterate")
save(dds, file='datasets/dds_iterated_size_factors.Rdata')
dds <- estimateDispersions(dds)
dds.t<-dds #for transplant
dds.o<-dds #for origin


#lrt test for transplant
dds.t<-nbinomLRT(dds.t, reduced=~origin)
resultsNames(dds.t)
traco=results(dds.t, contrast = c('transplant', 'O', 'K'), independentFiltering=F)
traco = traco[order(traco$pvalue),]
summary(traco)   #no differentially methylated genes for this test
head(na.omit(traco))
save(traco, file="datasets/traco_GENEBODIES_meth_p-1000_200.Rdata")

#lrt test for origin
dds.o<-nbinomLRT(dds.o, reduced=~transplant)
resultsNames(dds.o)
orico=results(dds.o, contrast = c('origin', 'O', 'K'), independentFiltering=F)
orico=orico[order(orico$pvalue),]
summary(orico)
head(na.omit(orico))
save(orico, file="datasets/orico_GENEBODIES_meth_p-1000_200.Rdata")


#save all the important objects
save(meth.rld, meth.coldata, dds.t, dds.o, counts, orico, traco,  file="datasets/deseqObjects_GENEBODIES_promoter1000_200.Rdata")




############ look at the distribution of mbd-scores for transplant and origin differential methylation ##########
top = 1000
YLIM=c(0, 0.25)
par(mfrow=c(1,2))
sig.o = na.omit(orico)[1:top,]
plot_subset_mbd_density(rownames(sig.o), 'red', YLIM = YLIM, MAIN = 'Density for Diff. Meth. by Origin')
sig.t = na.omit(traco)[1:top,]
plot_subset_mbd_density(rownames(sig.t), 'red', YLIM = YLIM, MAIN = 'Density for Diff. Meth. by Transplant')
###############################################################################################################


#---------------- save/load ----------------#
#load functions and data
directory<-"~/gitreps/reciprocal_transplant_methylation"
setwd(directory)
source("~/gitreps/reciprocal_transplant_methylation/scripts/reciprocal_methylation_project_functions.R")
lnames = load("datasets/ProteinTable.Rdata")
lnames = load("datasets/deseqObjects_GENEBODIES_promoter1000_200.Rdata")

#---------------- plot correlation heatmap ----------------#
library(pheatmap)
meth.rld.df = assay(meth.rld)
colnames(meth.rld.df) = colnames(counts)
labs = sub("_2m", "", colnames(meth.rld.df))
labs[labs=="KK4"]<-"KK4'"
labs[labs=="KO4"]<-"KO4'"
pheatmap(cor(meth.rld.df, method='pearson'), cluster_cols = T, cluster_rows = T, clustering_distance_rows = "maximum", clustering_distance_cols = "maximum", labels_row=labs, labels_col=labs)
dim(meth.rld.df)


############ OUTPUT FOR GO-MWU TESTS ############
lnames=load('datasets/ProteinTable.Rdata')
lnames
tot.uniq=length(unique(ptable$locusName))
ptable=ptable[!duplicated(ptable$locusName),]
tot.uniq
rownames(ptable) = ptable$locusName
nrow(ptable)
head(ptable)
o1=merge(data.frame(orico), ptable, by=0)
head(o1)
x<-o1$log2FoldChange < 0
x[x==T]<- -1
x[x==F]<- 1
stat=o1$stat*x
out = na.omit(data.frame(o1$genbank.prot, stat))
colnames(out)=c('prot', 'stat')
head(out)
write.csv(out, file='go_mwu/origin_gbm_go_mwu_input.csv', row.names=F, quote=F)



