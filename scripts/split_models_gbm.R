#split_models_gbm.R
#SET UP THE DATA TO RUN DESEQ
library('DESeq2')
#SAVE YOUR DIRECTORY NAME AND SET WD
directory<-"~/gitreps/reciprocal_transplant_methylation"
setwd(directory)


#load the data
lnames=load("datasets/gbm_counts_p-1000_200.Rdata")
head(gcounts)
dim(gcounts)

######## set up two different datasets, one with all timepoints and one with just 3-month samples (timepoint 2)
#gather data
all.samples=colnames(gcounts)
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
allcc = gcounts[,csample]
dim(allcc) #all samples, 44 3-month samples and 6 6-month samples
cc= allcc[,time=='2']
dim(cc)    #timepoint2 only samples
############################################################


#------------ set up DESeq2 inputs for both datasets and get 

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



#save/load
# save(counts,conditions,file="datasets/BayRT_remapped2genome_MBDseq.RData")
# lnames=load('BayRT_remapped2genome_MBDseq.RData')
# lnames


#get variance stabilized counts for timepoint 2
library('DESeq2')
dds<-DESeqDataSetFromMatrix(counts,
	colData = conditions, 
	design = formula(~ origin+transplant))
rl=rlog(dds)
vsd=assay(rl)
library(pheatmap)
pheatmap(cor(vsd))
save(dds,vsd,rl,conditions,counts,file="datasets/mbd_3mo.RData") 

#prepare conditions data for timepoint 2 and 3
head(allcc)
sample=colnames(allcc)

origin=substr(sample, start=1, stop=1)
transplant=substr(sample, start=2, stop=2)
num=substr(id[time=='2'], start=3, stop=5)
colony.id=paste(origin,num,sep="")
conditions=data.frame(sample,colony.id,origin,transplant, time)
means=apply(allcc,1,mean)
table(means>2)
allcounts=allcc[means>2,]


#get variance stabilized counts for timepoint 2 and 3
library('DESeq2')
dds<-DESeqDataSetFromMatrix(allcounts,
	colData = conditions, 
	design = formula(~ origin+transplant))
rl=rlog(dds)
vsd=assay(rl)
colnames(vsd) = colnames(allcounts)
library(pheatmap)
pheatmap(cor(vsd), method='maximum')
save(dds,vsd,rl,conditions,allcounts,file="mbd_3_and_6mo.RData") 






#now reload the 3month only data
rm(list=ls())
load("datasets/BayRT_remapped2genome_MBDseq.RData")

#----------------------
# splitting into subsets

lnames=load("datasets/mbd_3mo.RData")
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
o2kd=estimateSizeFactors(o2kd, type='iterate')
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
k2od=estimateSizeFactors(k2od, type='iterate')
k2od=estimateDispersions(k2od)
k2od=nbinomWaldTest(k2od)
resultsNames(k2od)
k2o.r=results(k2od, contrast = c('transplant', 'O', 'K'))
summary(k2o.r)

#origin effect for corals replanted at their original sites (KK and OO)
homed=estimateSizeFactors(homed, type='iterate')
homed=estimateDispersions(homed)
homed=nbinomWaldTest(homed)
resultsNames(homed)
home.r=results(homed, contrast = c('origin', 'O', 'K'))
summary(home.r)


#origin effect for corals placed at orpheus (KO and OO)
okatod=estimateSizeFactors(okatod, type='iterate')
okatod=estimateDispersions(okatod)
okatod=nbinomWaldTest(okatod)
resultsNames(okatod)
okato.r=results(okatod, contrast = c('origin', 'O', 'K'))
summary(okato.r)

#origin effect for corals placed at keppel (KK and OK)
okatkd=estimateSizeFactors(okatkd, type='iterate')
okatkd=estimateDispersions(okatkd)
okatkd=nbinomWaldTest(okatkd)
resultsNames(okatkd)
okatk.r=results(okatkd, contrast = c('origin', 'O', 'K'))
summary(okatk.r)


#origin effect for transplants
transd=estimateSizeFactors(transd, type='iterate')
transd=estimateDispersions(transd)
transd=nbinomWaldTest(transd)
resultsNames(transd)
trans.r=results(transd, contrast = c('origin', 'O', 'K'))
summary(trans.r)
save(trans.r, file='datasets/transr.Rdata')




save(okatk.r,okato.r,k2o.r,o2k.r,home.r,file="splitModels_MBD.RData")

#---------------------
