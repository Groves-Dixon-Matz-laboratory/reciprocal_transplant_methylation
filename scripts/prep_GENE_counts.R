#prep_GENE_counts.R


#---------------- Upload the data ------------------
#SET UP THE DATA TO RUN DESEQ
library('DESeq2')
#SAVE YOUR DIRECTORY NAME AND SET WD
directory<-"~/gitreps/reciprocal_transplant_methylation/"
setwd(directory)


counts=read.table('datasets/GENE_counts_3-6-18.tsv',header=TRUE,row.names=1,sep="\t"); head(counts)  #counts to full annotated gene regions
head(counts)
dim(counts)


#now remove the special cases counts
remove = rownames(counts)[grep("__", rownames(counts))]
counts = counts[!rownames(counts) %in% remove,]
dim(counts)


#remove other werid genes
keep = grep("LOC", rownames(counts))
length(keep)
rem=rownames(counts)[-keep]
print(rem)
counts=counts[keep,]


gcounts=counts
save(gcounts, file="datasets/gbm_counts_FULLGENE.Rdata")