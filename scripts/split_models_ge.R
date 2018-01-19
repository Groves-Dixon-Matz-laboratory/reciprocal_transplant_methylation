#split_models_ge.R
#SET UP THE DATA TO RUN DESEQ
library('DESeq2')
#SAVE YOUR DIRECTORY NAME AND SET WD
directory<-"~/gitreps/reciprocal_transplant_methylation"
setwd(directory)


#load the data
cc=read.table("datasets/all_rnaseq_counts.txt", header = T, row.names='geneID')
head(cc)
dim(cc)


#------------ set up DESeq2 inputs for both datasets and get 

sample=names(cc)
origin=substr(sample, start=1, stop=1)
transplant=substr(sample, start=2, stop=2)
num=substr(sample, start=3, stop=5)
colony.id=paste(origin,num,sep="")
conditions=data.frame(sample,colony.id,origin,transplant)
means=apply(cc,1,mean)
table(means>2)
counts=cc[means>2,]



#get variance stabilized counts for timepoint 2
library('DESeq2')
dds<-DESeqDataSetFromMatrix(counts,
	colData = conditions, 
	design = formula(~ origin+transplant))
rl=rlog(dds)
vsd=assay(rl)
library(pheatmap)
pheatmap(cor(vsd))
save(dds,vsd,rl,conditions,counts,file="datasets/ge_3mo_v2.RData") 


#----------------------
# splitting into subsets

lnames=load("datasets/ge_3mo_v2.RData")
library(DESeq2)
# library(empiricalFDR.DESeq2)
#sim=simulateCounts(dds)
#counts=assay(sim)
#head(counts)

home=which(conditions$origin==conditions$transplant) #subset of KK and OO samples
o2k=which(conditions$origin=="O")                    #subset for samples from Orpheus
k2o=which(conditions$origin=="K")                    #subset for samples from Keppel
okato=which(conditions$transplant=="O")              #subset for samples transplanted to Orpheus (OO and OK)
okatk=which(conditions$transplant=="K")              #subset for samples transplanted to Orpheus (KK and KO)

treat=paste(conditions$transplant, conditions$origin, sep='')
transplants=c(which(treat=="KO"), which(treat=="OK")) #subset for samples transplanted to Orpheus (KK and KO)


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

#test for origin effects in transplants (canalized transcription)
transd=DESeqDataSetFromMatrix(counts[,transplants],
	colData = conditions[transplants,], 
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


#origin effect for transplants
transd=estimateSizeFactors(transd)
transd=estimateDispersions(transd)
transd=nbinomWaldTest(transd)
resultsNames(transd)
trans.r=results(transd, contrast = c('origin', 'O', 'K'))
summary(trans.r)
save(trans.r, file='datasets/transr_GE.Rdata')




save(okatk.r,okato.r,k2o.r,o2k.r,home.r,file="datasets/splitModels_GE_v2.RData")

#---------------------
