#---------------- Upload the data ------------------
#SET UP THE DATA TO RUN DESEQ
library('DESeq2')
#SAVE YOUR DIRECTORY NAME AND SET WD
directory<-"~/gitreps/reciprocal_transplant_methylation/"
setwd(directory)

#READ IN THE COUNTS DATA. THESE ARE EXPORTED AT THE END OF THE MBD-seq data processing pipeline (see MBD-# seq_Data_Processing_Walkthrough).
counts=read.table('datasets/promoter_gbm_counts_p-1000_200.tsv',header=TRUE,row.names=1,sep="\t"); head(counts)  #small promo
head(counts)
dim(counts)


ubcounts = counts[,grep("ub", colnames(counts))]

counts = counts[,grep("2m", colnames(counts))]
head(counts)
dim(counts)



#get fraction for MBD
scount = counts[grep("__", rownames(counts)),]
gcounts = counts[-grep("__", rownames(counts)),]
head(gcounts)

genetots = apply(gcounts, 2, sum)
n.genetots = scount['__no_feature',]
amb.tots = scount['__too_low_aQual',] + scount['__ambiguous',]

sum(names(genetots) == names(n.genetots)) == length(genetots)
full.tot = genetots + n.genetots
frac.gene = genetots/full.tot
frac.not = n.genetots/full.tot
frac.amb = amb.tots/(full.tot + amb.tots)
mres = data.frame(as.numeric(frac.gene), 'MBD')
mamb = data.frame(as.numeric(frac.amb), 'MBD')
colnames(mres) = c('frac', 'type')
colnames(mamb) = c('frac', 'type')


grand.tot = full.tot + amb.tots
full.g.frac = genetots/grand.tot
mfres = data.frame(as.numeric(full.g.frac), 'MBD')
colnames(mfres) = c('frac', 'type')


#get fraction for UB
scount = ubcounts[grep("__", rownames(ubcounts)),]
gcounts = ubcounts[-grep("__", rownames(ubcounts)),]
amb.tots = scount['__too_low_aQual',] + scount['__ambiguous',]
head(gcounts)

genetots = apply(gcounts, 2, sum)
n.genetots = scount['__no_feature',]
sum(names(genetots) == names(n.genetots)) == length(genetots)
full.tot = genetots + n.genetots

frac.gene = genetots/full.tot
frac.not = n.genetots/full.tot
frac.amb = amb.tots/(full.tot + amb.tots)
ures = data.frame(as.numeric(frac.gene), 'flowthrough')
uamb = data.frame(as.numeric(frac.amb), 'flowthrough')
colnames(ures) = c('frac', 'type')
colnames(uamb) = c('frac', 'type')


grand.tot = full.tot + amb.tots
full.g.frac = genetots/grand.tot
ufres = data.frame(as.numeric(full.g.frac), 'flowthrough')
colnames(ufres) = c('frac', 'type')


res=rbind(mres, ures)
amb=rbind(mamb, uamb)
fres = rbind(mfres, ufres)






boxplot(res[,1]*100~res[,2], ylab="% of Mapped Reads mapping to annotated CDS")
boxplot(amb[,1]*100~amb[,2], ylab="% of Abmiguously Mapped Reads")
boxplot(fres[,1]*100~fres[,2], ylab="% of All Reads mapping to annotated CDS")

points(res[,1]*100 ~ jitter(rep(1, nrow(gres)), factor=2))
nres = data.frame(as.numeric(frac.not), 'mapped to coding')





