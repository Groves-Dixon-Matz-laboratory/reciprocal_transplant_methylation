#doublecheck_dapc.R
#this script does the same things as DAPC_gbm_transplant.R
#recoding it to make sure there are no bugs.



#### SETUP ####
library(adegenet)
setwd("~/gitreps/reciprocal_transplant_methylation/")
source("~/gitreps/reciprocal_transplant_methylation/scripts/reciprocal_methylation_project_functions.R")



##### UPLOAD DATA ####
ll=load("datasets/splitModels_MBD.RData")
ll=load("datasets/mbd_3mo.RData")
head(k2o.r)
head(o2k.r)

######## BUILD VOLCANO PLOTS ########
#these illustrate differential methylation within populations due to transplantation
par(mfrow=c(1,1))
YLIM=c(0,16)
volcano_plot(k2o.r, XLIM=c(-1.5, 1.5), YLIM=YLIM, MAIN='KK vs KO', draw.box=F) #KK vs KO test
volcano_plot(o2k.r, XLIM=c(-1.5, 1.5), YLIM=YLIM, MAIN='OO vs OK', draw.box=F) #OO vs OK test
par(mfrow=c(1,1))



#assign plastic genes
CUT=0.01
s1=k2o.r[k2o.r$pvalue<CUT,]
nrow(s1)
s2=o2k.r[o2k.r$pvalue<CUT,]
nrow(s2)
sg1=rownames(s1)
sg2=rownames(s2)
sg=unique(c(sg1, sg2))
length(sg)

vsds=vsd[sg,]
nrow(vsds)
#560

#fix samples names
s=colnames(vsds)
fs=sub("_2m", "", s)
colnames(vsds) = fs

######## SET SPLIT DATASTS FOR NATIVE CORALS AND TRANSPLANTED CORALS ########
samples = colnames(vsds)
natives = append(samples[grep("KK", samples)], samples[grep("OO", samples)])
transplants = append(samples[grep("KO", samples)], samples[grep("OK", samples)])
n.vsd = vsds[,natives]
t.vsd = vsds[,transplants]


######## PREPARE DATA FOR DAPC ########
tn.vsd=t(n.vsd)
pcp=prcomp(tn.vsd, retx=T, center=T, scale=T)
scores=pcp$x
clus=find.clusters(tn.vsd, max.n.clus=15, n.clust=2, n.pca=4)
n.treats = substr(colnames(n.vsd), start=1, stop=2)
clus$grp=n.treats

######## BUILD DISCIMINATE FUNCTION ########
dp=dapc(tn.vsd,clus$grp, n.da=1, perc.pca=80) 


#plot distributions
quartz()
scatter(dp,bg="white",scree.da=FALSE,legend=TRUE,solid=.4)


######## APPLY DISCIMINATE FUNCTION TO TRANSPLANTED SAMPLES ########
tt.vsd = t(t.vsd)
tdp = predict.dapc(dp, newdata=tt.vsd)
names(tdp)
tdp$ind.scores

dim(vsds)
all.dp = predict.dapc(dp, newdata=t(vsds))
a=all.dp$ind.scores



#CHECK THAT EVERYTHING MATCHES
#for transplants
rownames(tdp$ind.scores)
tt = a[rownames(tdp$ind.scores),]
length(tt)
tt == tdp$ind.scores
#for natives
names(dp)
x=dp$ind.coord
rownames(x)
tn=a[rownames(x),]
tn
x
names(tn)==rownames(x)
y=data.frame(tn)
x=data.frame(x)
colnames(y) = c('LD1')
x==y
z=x$LD1 - y$LD1
z


######## CONNECT DAPC RESULTS WITH PHENOTYPIC DATA ########

#UPLOAD THE TRAIT DATA TO MERGE WITH MBD-SEQ RESULTS
lnames=load('datasets/bayRT_rlog_conditions_traits.RData')
traits
a=data.frame(a)
a$Colony.ID = rownames(a)
a

x=merge(traits, a, by = 'Colony.ID', all.x=T)
dim(x)
dim(traits)
head(x)
sum(x$Colony.ID==traits$Colony.ID) == nrow(traits) #check it's lined up
traits$gbm.ld1 = x$LD1
head(traits)


#Weird lipid measure
hist(traits$LIPID)

a$treat = substr(a$Colony.ID, start = 1, stop=2)
head(a)

get.colors = function(treat){
	c=treat
	c[c=="KK"]<-'blue'
	c[c=="KO"]<-'cyan'
	c[c=="OK"]<-'orange'
	c[c=="OO"]<-'red'
	return(c)
}


#look at correlations
head(traits)
traits$treat = paste(traits$ori, traits$tra, sep='')
traits$colors=get.colors(traits$treat)
head(traits)

fp="GAIN"
plot(traits[,fp]~traits[,'gbm.ld1'], pch=21, bg=traits$colors, cex=2)
for (t in unique(traits$treat)){
	sub=traits[traits$treat == t,]
	lms=lm(sub[,fp]~ sub[,'gbm.ld1'])
	print("----------------")
	print(paste("Results for", t))
	print(summary(lms))
	abline(lms, col=sub$colors[1])
}

#get convergence scores
ntraits = traits[traits$Colony.ID %in% natives,]
nmns = tapply(ntraits$gbm.ld1, INDEX=ntraits$treat, mean)
abline(v= nmns)
ttraits=traits[traits$Colony.ID %in% transplants,]
traits$target = NA
traits$target[traits$tra == 'O']<-nmns['OO']
traits$target[traits$tra == 'K']<-nmns['KK']

#subset for transplants
ttraits=traits[traits$Colony.ID %in% transplants,]
ttraits$distance=abs(ttraits$target - ttraits$gbm.ld1)
mnd = mean(ttraits$distance)
sdd = sd(ttraits$distance)
ttraits$conv = (mnd - ttraits$distance)/sdd
ttraits


#plot relationships
fp="GAIN"
plot(ttraits[,fp]~ttraits[,'conv'], pch=21, bg=ttraits[,'colors'], cex=2)
lmf = lm(ttraits[,fp]~ttraits[,'conv'])
summary(lmf)

par(mfrow = c(3,2))
for (fp in c('GAIN', 'LIPID', 'CARB', 'PROTEIN', 'ZOOX')){
	lmf = lm(ttraits[,fp]~ttraits[,'conv'])
	p=round(summary(lmf)$coefficients[2,4], digits=2)
	r=round(summary(lmf)$r.squared, digits=2)
	plot(ttraits[,fp]~ttraits[,'conv'], pch=21, bg=ttraits[,'colors'], cex=2, ylab=fp, xlab='convergence')
	title(main=bquote("R2="*.(r)))
	title(main=bquote('p='*.(p)), line=1)
	abline(lmf, lwd=2, lty=2)
	print("-------------------")
	print(paste("Results for", fp))
	print(summary(lmf))
}






