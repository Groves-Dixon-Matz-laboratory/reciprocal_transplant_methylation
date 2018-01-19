#---------------- Upload the data ------------------
#SET UP THE DATA TO RUN DESEQ
library('DESeq2')
#SAVE YOUR DIRECTORY NAME AND SET WD
directory<-"~/gitreps/reciprocal_transplant_methylation"
setwd(directory)

cc=read.table("all_rnaseq_counts.txt",header=T,sep="\t")
head(cc)
row.names(cc)=cc$geneID
cc=cc[,-1]

sample=names(cc)
origin=sub("[OK][0-9]+","",sample)
transplant=sub("^[OK]","",sample)
transplant=sub("[0-9]+","",transplant)
num=sub("[OK]+","",sample)
colony.id=paste(origin,num,sep="")
conditions=data.frame(cbind(sample,colony.id,origin,transplant))
means=apply(cc,1,mean)
table(means>2)
counts=cc[means>2,]
save(counts,conditions,file="BayRT_remapped2genome_tagseq.RData")
load("BayRT_remapped2genome_tagseq.RData")

library('DESeq2')
dds<-DESeqDataSetFromMatrix(counts,
	colData = conditions, 
	design = formula(~ origin+transplant))
rl=rlog(dds)
vsd=assay(rl)
library(pheatmap)
pheatmap(cor(vsd))
save(dds,vsd,rl,conditions,counts,file="datasets/GE_3mo.RData") 

#----------------------
# splitting into subsets

load("datasets/GE_3mo.RData")
library(DESeq2)
# library(empiricalFDR.DESeq2)
#sim=simulateCounts(dds)
#counts=assay(sim)
#head(counts)

home=which(conditions$origin==conditions$transplant)
o2k=which(conditions$origin=="O")
k2o=which(conditions$origin=="K")
okato=which(conditions$transplant=="O")
okatk=which(conditions$transplant=="K")

homed=DESeqDataSetFromMatrix(counts[,home],
	colData = conditions[home,], 
	design = formula(~ origin))

o2kd=DESeqDataSetFromMatrix(counts[,o2k],
	colData = conditions[o2k,], 
	design = formula(~ colony.id+transplant))

k2od=DESeqDataSetFromMatrix(counts[,k2o],
	colData = conditions[k2o,], 
	design = formula(~ colony.id+transplant))

okatod=DESeqDataSetFromMatrix(counts[,okato],
	colData = conditions[okato,], 
	design = formula(~ origin))

okatkd=DESeqDataSetFromMatrix(counts[,okatk],
	colData = conditions[okatk,], 
	design = formula(~ origin))

o2kd=DESeq(o2kd)
o2k.r=results(o2kd)
summary(o2k.r)

k2od=DESeq(k2od)
k2o.r=results(k2od)
summary(k2o.r)

homed=DESeq(homed)
home.r=results(homed)
summary(home.r)

okatod=DESeq(okatod)
okato.r=results(okatod)
summary(okato.r)

okatkd=DESeq(okatkd)
okatk.r=results(okatkd)
summary(okatk.r)

save(okatk.r,okato.r,k2o.r,o2k.r,home.r,file="splitModels_GE.RData")

#--------------------- start here to do transplant figures ------------
directory<-"~/gitreps/reciprocal_transplant_methylation"
setwd(directory)
ll=load("datasets/splitModels_GE.RData")
ll=load("datasets/GE_3mo.RData")
gehome=home.r
geo2k=o2k.r
gek2o=k2o.r
geokatk=okatk.r
geokato=okato.r
gevsd=vsd
geconditions=conditions


ll=load("datasets/splitModels_MBD.RData")
ll=load("datasets/mbd_3mo.RData")
conditions$sample=sub("_2m","",conditions$sample)
# sshared=intersect(geconsditions$sample,conditions$sample)
vsd=data.frame(vsd)
colnames(vsd) = conditions$sample
gevsd=data.frame(gevsd)
colnames(gevsd) = geconditions$sample



library(vegan)
ad=adonis(t(gevsd)~origin+colony.id+transplant,data= geconditions,method="manhattan")
labs=c("origin","colony.id","transplant","residuals")
cols=c("skyblue","green2","coral","grey80")
cols = append(gg_color_hue(3), 'grey')
labs2 = paste(labs, round(ad$aov.tab$R2[1:4]*100, digits=1))
pie(ad$aov.tab$R2[1:4],labels=labs2,col=cols,main="Gene expression")

ad=adonis(t(vsd)~origin+colony.id+transplant,data= conditions,method="manhattan")
labs2 = paste(labs, round(ad$aov.tab$R2[1:4]*100, digits=1))
pie(ad$aov.tab$R2[1:4],labels=labs2,col=cols,main="Methylation")

#------------------

########## use this chunk for gene expression #########
# vsds=gevsd[abs(gehome$stat)>2.5,]
vsds=gevsd[abs(gehome$pvalue)<0.01,]
coo=geconditions
colnames(vsds) = coo$sample
str(gehome)
dim(vsds)
#######################################################

#### use this chunk for 3 month only methylation ######
vsds=vsd[abs(home.r$stat)>2.5,]
vsds=vsd[!is.na(home.r$pvalue) & abs(home.r$pvalue)<0.01,]
vsds=vsd[abs(home.r$log2FoldChange) > log(1.5, 2),]
coo=conditions
str(vsds)
dim(vsds)
# GE: 1186 genes
# MBD: 398 genes
#######################################################

library(vegan)
#library(rgl)
library(ape)
library(ggplot2)
library(MASS)

#set up PCoA
dd.pcoa=pcoa(vegdist(t(vsds),method="manhattan"))
scores=dd.pcoa$vectors
moo=apply(scores[grep("OO",coo$sample),],2,mean)
mok=apply(scores[grep("OK",coo$sample),],2,mean)
mkk=apply(scores[grep("KK",coo$sample),],2,mean)
mko=apply(scores[grep("KO",coo$sample),],2,mean)

# plotting PCoA
margin=5
pc1=1;pc2=2
plot(scores[, pc1], scores[, pc2],xlab="PCo1",ylab="PCo2",mgp=c(2.3,1,0),xlim=c(min(scores[,1])-margin,max(scores[,1]+margin)),ylim=c(min(scores[,2])-margin,max(scores[,2]+margin)),cex=0.3,pch=19)
#ordispider(scores[, c(pc1, pc2)],conditions$colony.id,col="skyblue")
o2o=grep("OO", coo$sample)
o2k=grep("OK", coo$sample)
arrows(scores[o2o , pc1],scores[o2o, pc2] ,scores[o2k , pc1],scores[o2k, pc2],length=0.07)
k2o=grep("KO", coo$sample)
k2k=grep("KK", coo$sample)
arrows(scores[k2k , pc1],scores[k2k, pc2] ,scores[k2o , pc1],scores[k2o, pc2],length=0.07)
#ordihull(scores[, c(pc1, pc2)], coo$origin,draw="polygon",col="grey90",label=T,cex=0.8)
#ordiellipse(scores[, c(pc1, pc2)], coo$origin,label=T,draw="polygon",col="grey90",cex=1)
ordihull(scores[coo$origin=="K", c(pc1, pc2)],draw="polygon",col="skyblue",label=T,cex=0.8,groups=coo$origin)
ordihull(scores[coo$origin=="O", c(pc1, pc2)],draw="polygon",col="pink",label=T,cex=0.8,groups=coo$origin)
arrows(moo[pc1],moo[pc2],mok[pc1],mok[pc2],length=0.15,lwd=3,col="red")
arrows(mkk[pc1],mkk[pc2],mko[pc1],mko[pc2],length=0.15,lwd=3,col="blue")




# lnames=load('/Users/grovesdixon/git_Repositories/reciprocal_transplant_methylation/datasets/deseqObjects_GENEBODIES_promoter1000_200.Rdata')
# lnames
# vsd=data.frame(assay(meth.rld))
# colnames(vsd) = sub("_2m", "", colnames(counts))
# head(vsd)
# dim(vsd)
# head(traco)
# sum(rownames(traco)==rownames(vsd))==nrow(vsd) #check rownames are same
# vsds=vsd[abs(traco$stat>2.5),]
# vsds=na.omit(vsds)
# dim(vsds)



################### use this chunk for methylation ################
head(k2o.r)
head(o2k.r)
rownames(k2o.r) = rownames(home.r)
rownames(o2k.r) = rownames(home.r)
CUT=2.5
# top=1000
# #for transplant among keppel corals
# sig.genes = rownames(k2o.r)[abs(k2o.r$stat)>2.2];length(sig.genes)
# sig.genes = tail(rownames(k2o.r[order(abs(k2o.r$stat)),]), n=top);length(sig.genes)

# #for transplant among orpheus corals
# sig.genes = rownames(o2k.r)[abs(o2k.r$stat)>3.2];length(sig.genes)
# sig.genes = tail(rownames(o2k.r[order(abs(o2k.r$stat)),]), n=top);length(sig.genes)


#for both
sig.genes1 = rownames(k2o.r)[abs(k2o.r$stat)>CUT];length(sig.genes1)
sig.genes2 = rownames(o2k.r)[abs(o2k.r$stat)>CUT];length(sig.genes2)
sig.genes = append(sig.genes1, sig.genes2)

#to use top set significant genes:
# sig.genes1 = tail(rownames(k2o.r[order(abs(k2o.r$stat)),]), n=top);length(sig.genes)
# sig.genes2 = tail(rownames(o2k.r[order(abs(o2k.r$stat)),]), n=top);length(sig.genes)
# sig.genes = append(sig.genes1, sig.genes2)
# length(sig.genes)


sig.genes=unique(sig.genes)
sig.genes=sig.genes[!is.na(sig.genes)]
head(sig.genes)
length(sig.genes)
vsds=vsd[sig.genes,]
coo=conditions
colnames(vsds)=sub("_2m", "", colnames(vsds))
###########################output for GO
not.sig = rownames(k2o.r)[!rownames(k2o.r) %in% sig.genes]
length(not.sig)+length(sig.genes) == nrow(k2o.r)
ns = data.frame(not.sig, 0)
sig=data.frame(sig.genes, 1)
colnames(ns) = c('locusName', 'sig')
colnames(sig) = c('locusName', 'sig')
out=rbind(sig, ns)
head(out)
lnames=load("~/git_Repositories/reciprocal_transplant_methylation/datasets/ProteinTable.Rdata")
head(ptable)
p=merge(ptable, out, by='locusName')
head(p)
out=p[,c('genbank.prot', 'sig')]
head(out)
write.csv(out, file='~/git_Repositories/reciprocal_transplant_methylation/go_mwu/plastic_meth_genes.csv', row.names=F, quote=F)


############################where are they in the MBD-score distribution?
lnames=load("~/git_Repositories/reciprocal_transplant_methylation/datasets/capturedVflowthroughResults_p-1000_200.Rdata")
head(res)
hist(res$log2FoldChange, breaks=30)
sub=res[sig.genes,]
dim(sub)
hist(sub$log2FoldChange, add=T, col='blue', breaks=30)
plot(density(sub$log2FoldChange), col='red')
lines(density(na.omit(res$log2FoldChange)))

###################################################################

################### use this chunk for gene expression ################
head(gek2o)
head(geo2k)
sig.genes = rownames(gek2o)[abs(gek2o$stat)>4]
sig.genes = append(sig.genes, rownames(geo2k)[abs(geo2k$stat)>2.5])
sig.genes=unique(sig.genes)
sig.genes=sig.genes[!is.na(sig.genes)]
head(sig.genes)
length(sig.genes)
vsds=gevsd[sig.genes,]
coo=geconditions
###################################################################



##set of discriminant analysis
lin = data.frame(t(vsds))
dim(lin)
ori=substr(rownames(lin), start=1, stop=1)
trans=substr(rownames(lin), start=2, stop=2)
treat=substr(rownames(lin), start=1, stop=2)
# hometags=c('OO', 'KK')
# home=treat
# home[!home %in% hometags]<-NA



# lin$ori=ori
# lda <- lda(ori ~ ., data=lin)


lin$trans=trans
lda <- lda(trans ~ ., data=lin)


lda.vals <- predict(lda)
# ldahist(data = lda.vals$x[,1], g=treat)
lda.res = data.frame(lda.vals$x, treat)
num=substr(rownames(lda.res), start=3, stop=5)
lda.res=lda.res[!num=='10',]

mns=tapply(lda.res$LD1, lda.res$treat, mean)

#plot discriminant analysis
color.set=c('blue', 'cyan', 'orange', 'red')
quartz()
ggplot(data=lda.res, aes(LD1, fill=treat, color=treat)) + geom_density(alpha=0.4) + theme(panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +scale_color_manual(values=color.set) + scale_fill_manual(values=color.set) # + geom_vline(xintercept=mns, color=color.set)


tlda=lda.res

head(tlda)
head(lda.res)
lda.res$LDt = tlda$LD1

plot(LD1~LDt, data=lda.res, xlab='LD transplant', ylab='LD origin')
points(LD1~LDt, data=lda.res[lda.res$treat=='KK',], bg='cyan', pch=21, cex=2)
points(LD1~LDt, data=lda.res[lda.res$treat=='KO',], bg='blue', pch=21, cex=2)
points(LD1~LDt, data=lda.res[lda.res$treat=='OK',], bg='orange', pch=21, cex=2)
points(LD1~LDt, data=lda.res[lda.res$treat=='OO',], bg='red', pch=21, cex=2)
legend('topleft', c('KK', 'KO', 'OK', 'OO'), pt.bg=c('cyan', 'blue', 'orange', 'red'), pch=21)


#============
# addiing 1st PCoAs to traits

lnames=load('bayRT_rlog_conditions_traits.RData')
lnames=load('~/Documents/GRANTS/IOS2016/genotype_reps_included/mcmc_MBD_3mo/splitModels_GE_MBD.RData')

vsds=vsd[abs(home.r$stat)>2.5,]
library(vegan)
#library(rgl)
library(ape)
dd.pcoa=pcoa(vegdist(t(vsds),method="manhattan"))
scores=dd.pcoa$vectors
row.names(scores)=sub("_2m","",row.names(scores))

#gather the mbd.pco1 data and add to the conditions
mbd.pco1=c()
for(co in as.character(traits$Colony.ID)){
	if (co %in% row.names(scores)) {
		mbd.pco1=append(mbd.pco1,scores[co,1])
	} else { mbd.pco1=append(mbd.pco1,NA) }
}
traits$mbd.pco1=mbd.pco1
#gather the gene expression pco1 data and add to the conditions
vsds=gevsd[abs(gehome$stat)>2.5,]
dd.pcoa=pcoa(vegdist(t(vsds),method="manhattan"))
scores=dd.pcoa$vectors
ge.pco1=c()
for(co in as.character(traits$Colony.ID)){
	if (co %in% row.names(scores)) {
		ge.pco1=append(ge.pco1,scores[co,1])
	} else { ge.pco1=append(ge.pco1,NA) }
}
traits$ge.pco1=ge.pco1

quartz()
#--------gather the mbd discim data ------------
lda.res$Colony.ID = rownames(lda.res)
head(lda.res)
x = merge(traits, lda.res, by='Colony.ID', all.x=T)
x=x[order(x$Colony.ID),]
sum(x$Colony.ID == traits$Colony.ID)==nrow(traits)
traits$mbd.LD1=x$LD1
#------------------------------------------------


head(traits)
#plot correlation with gain
CEX=2
LWD=2
plot(GAIN~mbd.LD1, data=traits)
#KO samples
d=traits[traits$ori=='K' & traits$tra=='O' & !is.na(traits$mbd.LD1) & !is.na(traits$GAIN),]
points(GAIN~mbd.LD1, data=d, bg='cyan', pch=21, cex=CEX) #KO samples
lmko=lm(GAIN~mbd.LD1, data=d)
clip(x1=min(na.omit(d$mbd.LD1)), max(d$mbd.LD1), 0, 100)
abline(lmko, col = 'cyan', lty=2, lwd=LWD)
summary(lmko)
clip(-1e4,1e4,-1e4,1e4)
#OK samples
d=traits[traits$ori=='O' & traits$tra=='K' & !is.na(traits$mbd.LD1) & !is.na(traits$GAIN),]
points(GAIN~mbd.LD1, data=d, bg='orange', pch=21, cex=CEX)  #OK samples
lmok=lm(GAIN~mbd.LD1, data=d)
clip(x1=min(na.omit(d$mbd.LD1)), max(d$mbd.LD1), 0, 100)
abline(lmok, col = 'orange', lty=2, lwd=LWD)
summary(lmok)
clip(-1e4,1e4,-1e4,1e4)
#OO samples
d=traits[traits$ori=='O' & traits$tra=='O' & !is.na(traits$mbd.LD1) & !is.na(traits$GAIN),]
points(GAIN~mbd.LD1, data=d, bg='red', pch=21, cex=CEX)  #OO samples
lmoo=lm(GAIN~mbd.LD1, data=d)
clip(x1=min(na.omit(d$mbd.LD1)), max(d$mbd.LD1), 0, 100)
abline(lmoo, col = 'red', lty=2, lwd=LWD)
summary(lmoo)
clip(-1e4,1e4,-1e4,1e4)
#KK samples
d=traits[traits$ori=='K' & traits$tra=='K' & !is.na(traits$mbd.LD1) & !is.na(traits$GAIN),]
points(GAIN~mbd.LD1, data=d, bg='blue', pch=21, cex=CEX)
lmkk=lm(GAIN~mbd.LD1, data=d)
clip(x1=min(na.omit(d$mbd.LD1)), max(d$mbd.LD1), 0, 100)
abline(lmkk, col = 'blue', lty=2, lwd=LWD)
summary(lmkk)
clip(-1e4,1e4,-1e4,1e4)
legend('topright', c('KK', 'KO', 'OK', 'OO'), pt.bg=c('blue', 'cyan', 'orange', 'red'), pch=21)
head(traits)



#spider the genotypes
traits$geno=paste(traits$ori, traits$num, sep='')
x=na.omit(traits[,c('geno', 'mbd.LD1', 'GAIN')])
ordispider(x[, c(2,3)],x$geno,col="grey", lwd=.75)



summary(lmkk)
summary(lmko)
summary(lmok)
summary(lmoo)

#overlay mortality
#get mortality data
lnames=load('~/git_Repositories/reciprocal_transplant_methylation/datasets/wgcnaInitialize_p1000_200_iterate_mnCount10.Rdata')
head(traitData)
z=traitData[-grep('3m', rownames(traitData)),]
x=data.frame(z$Colony.ID, z$mortality)
colnames(x) = c('Colony.ID', 'mortality')
head(x)
m=merge(traits, x, by = 'Colony.ID', all.x=T)
sum(m$Colony.ID==traits$Colony.ID)==nrow(traits)
traits$mortality = m$mortality
#plot death points
head(traits)
dead=traits[traits$mortality==1 & !is.na(traits$mortality),]
points(GAIN~mbd.LD1, data=dead, pch=4, cex=2.5)
traits$geno=paste(traits$ori, traits$num, sep='')
dead=traits[traits$mortality==1 & !is.na(traits$mortality),]
dead.genos=dead$geno
died=traits[traits$geno %in% dead.genos,]

points(GAIN~mbd.LD1, data=died, pch=4, cex=2.5)
head(traits)
m=traits[traits$mortality==1 & !is.na(traits$mortality),]
m
mcols=paste(m$ori, m$tra, sep="")
mcols[mcols=='KK']<-'blue'
mcols[mcols=='KO']<-'cyan'
mcols[mcols=='OK']<-'orange'
mcols[mcols=='OO']<-'red'

abline(v=m$mbd.LD1, col=mcols)




#plot correlation for ge with gain
CEX=2
plot(GAIN~ge.LD1, data=traits)
points(GAIN~ge.LD1, data=traits[traits$ori=='K' & traits$tra=='O',], bg='blue', pch=21, cex=CEX) #KO samples
lmko=lm(GAIN~ge.LD1, data=traits[traits$ori=='K' & traits$tra=='O',]); abline(lmko, col = 'blue', lty=2);summary(lmko)
points(GAIN~ge.LD1, data=traits[traits$ori=='O' & traits$tra=='K',], bg='orange', pch=21, cex=CEX)  #OK samples
lmok=lm(GAIN~ge.LD1, data=traits[traits$ori=='O' & traits$tra=='K',]); abline(lmok, col = 'orange', lty=2);summary(lmok)
points(GAIN~ge.LD1, data=traits[traits$ori=='O' & traits$tra=='O',], bg='red', pch=21, cex=CEX)  #OO samples
lmoo=lm(GAIN~ge.LD1, data=traits[traits$ori=='O' & traits$tra=='O',]); abline(lmoo, col = 'red', lty=2);summary(lmoo)
points(GAIN~ge.LD1, data=traits[traits$ori=='K' & traits$tra=='K',], bg='cyan', pch=21, cex=CEX)  #KK samples
lmkk=lm(GAIN~ge.LD1, data=traits[traits$ori=='K' & traits$tra=='K',]); abline(lmkk, col = 'cyan', lty=2);summary(lmkk)
legend('topright', c('KK', 'KO', 'OK', 'OO'), pt.bg=c('cyan', 'blue', 'orange', 'red'), pch=21)


save("expr","conditions","traits",file="~/Documents/LineBay_RT_reanalysis/bayRT_rlog_conditions_traits.RData")
load("~/Documents/LineBay_RT_reanalysis/bayRT_rlog_conditions_traits.RData")

plot(mbd.pco1~ge.pco1,traits,ylab="Methylation PCo1", xlab="Transcription PCo1",mgp=c(2.1,1,0))
summary(lm(mbd.pco1~ge.pco1,traits))
abline(lm(mbd.pco1~ge.pco1,traits),col="red")

plot(GAIN~mbd.pco1,traits,ylab="daily weight gain, %", xlab="Methylation PCo1",mgp=c(2.1,1,0),col="white")
points(GAIN~mbd.pco1,traits[grep("KO",traits$Colony.ID),],col="blue",pch=19)
points(GAIN~mbd.pco1,traits[grep("KK",traits$Colony.ID),],pch=1,col="grey50")
clip(-100,0,0.05,0.5)
abline(lm(GAIN~mbd.pco1,traits[grep("KO",traits$Colony.ID),]),col="blue")
clip(-100,120,0,0.5)
summary(lm(GAIN~mbd.pco1,traits[grep("KO",traits$Colony.ID),])) # p = 0.25, R2=0.067
points(GAIN~mbd.pco1,traits[grep("OO",traits$Colony.ID),],col="grey80",pch=19)
points(GAIN~mbd.pco1,traits[grep("OK",traits$Colony.ID),],col="red",pch=19)
clip(-15,90,0,0.5)
abline(lm(GAIN~mbd.pco1,traits[grep("OK",traits$Colony.ID),]),col="red")
summary(lm(GAIN~mbd.pco1,traits[grep("OK",traits$Colony.ID),])) # p = 0.116, R2=0.19
clip(-10,120,0,0.5)
legend("topright",pch=c(19,19,1,19),col=c("blue","red","grey50","grey80"),c("K2O","O2K","K2K","O2O"),cex=0.8)

ttt=traits[grep("KO|OK",traits$Colony.ID),]
li=lm(GAIN~ori+abs(mbd.pco1),ttt)
summary(li)
plot(GAIN~abs(mbd.pco1),ttt)

mbdDist=abs(ttt$mbd.pco1[grep("KO",ttt$Colony.ID)]-moo[1])
mbdDist=append(mbdDist,abs(ttt$mbd.pco1[grep("OK",ttt$Colony.ID)]-mkk[1]))
ttt$mbdDist=mbdDist
plot(GAIN~mbdDist,ttt)
li=lm(GAIN~ori+mbdDist,ttt)
summary(li)
              # Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.2637653  0.0618251   4.266 0.000591 ***
# oriO        -0.0398432  0.0292927  -1.360 0.192631    
# mbdDist     -0.0010548  0.0005542  -1.903 0.075157 .  

sum(!is.na(traits$ge.pco1))

plot(GAIN~ge.pco1,traits)
points(GAIN~ge.pco1,traits[grep("KO",traits$Colony.ID),],col="skyblue",pch=19)
points(GAIN~ge.pco1,traits[grep("OK",traits$Colony.ID),],col="coral",pch=19)


#------------------------ descriminant analysis -----------------------


