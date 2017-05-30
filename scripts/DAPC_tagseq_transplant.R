#DAPC_tagseq_transplant.R
#run discriminate analysis of principal components on tag-seq data
#Groves Dixon
#last updated 5-1-17

#### SETUP ####
library(adegenet)
source("~/gitreps/reciprocal_transplant_methylation/scripts/reciprocal_methylation_project_functions.R")



##### UPLOAD DATA ####
setwd("~/gitreps/reciprocal_transplant_methylation/")
ll=load("datasets/splitModels_GE.RData") 
ll=load("datasets/GE_3mo.RData")
colnames(vsd) = colnames(counts)
head(k2o.r)
head(o2k.r)


######## BUILD VOLCANO PLOTS ########
#these illustrate differential expression within populations due to transplantation
par(mfrow=c(1,2))
YLIM=c(0, 12)
XLIM=c(-1.5, 1.5)
plot.k2o=data.frame(k2o.r)
plot.k2o= plot.k2o[plot.k2o $log2FoldChange<XLIM[2],]
volcano_plot(plot.k2o, XLIM=c(-1.5, 1.5), YLIM=YLIM, MAIN='KK vs KO', draw.box=F) #KK vs KO test
volcano_plot(o2k.r, XLIM=c(-1.5, 1.5), YLIM=YLIM, MAIN='OO vs OK', draw.box=F) #OO vs OK test
par(mfrow=c(1,1))


######### SUBSET FOR GENES THAT SHOW SIGNIFICANT DIFFERENCES DUE TO TRANSPLANTATION ########
#get count for FDR passing genes
CUT=.1
sig.genes1 = rownames(k2o.r)[!is.na(k2o.r$padj) & abs(k2o.r$padj)<CUT];length(sig.genes1) #sig for KK vs KO
sig.genes2 = rownames(o2k.r)[!is.na(o2k.r$padj) & abs(o2k.r$pvalue)<CUT];length(sig.genes2) #sig for OO vs OK
sig.genes = append(sig.genes1, sig.genes2)
sig.genes=unique(sig.genes)
sig.genes=sig.genes[!is.na(sig.genes)]
head(sig.genes)
length(sig.genes)


#actually select the genes
CUT=0.01  #unadjusted p-value cutoff
sig.genes1 = rownames(k2o.r)[abs(k2o.r$pvalue)<CUT];length(sig.genes1) #sig for KK vs KO
sig.genes2 = rownames(o2k.r)[abs(o2k.r$pvalue)<CUT];length(sig.genes2) #sig for OO vs OK
sig.genes = append(sig.genes1, sig.genes2)
sig.genes=unique(sig.genes)
sig.genes=sig.genes[!is.na(sig.genes)]
head(sig.genes)
length(sig.genes)
vsds=vsd[sig.genes,]
coo=conditions







#optionall set data for all genes
vsds=vsd


#look at subsetted dataset of variance stabilized counts
dat=vsds
degs10<-rownames(dat)
head(dat)
dim(dat) #1326 genes most differentially methylated by origin

#split variance stabilized counts into 'home' corals (KK and OO) and transplants (KO and OK)
a.vsd<-dat[,c(grep("OO",colnames(dat)), grep("KK",colnames(dat)))]     #home corals
a.vsd.supp<-dat[,c(grep("OK",colnames(dat)),grep("KO",colnames(dat)))] #transplants


#######################Some genes have insufficient variance for DFA in data subsets!...find those genes using loop below
dframe=(a.vsd[degs10,])
for(col in rownames(dframe))
	{  min=min(dframe[col,])
		max=max(dframe[col,])
		if(min == max){print(col)}
		}

#If print out above, copy gene names below and run line to remove genes from analysis. 

######################now use appropriate dataset for analysis

pcp=prcomp(t(a.vsd[degs10,]), retx=TRUE, center=TRUE, scale.=TRUE) #scale.=TRUE
scores=pcp$x
screeplot(pcp,bstick=T) # only 1st PC is non-random when using diff exp genes only; higher PCs non-random when using all data...

# adegenet: finding clusters (even though we know what clusters we want) - choose 4 PCs and 2 groups
clus=find.clusters(t(a.vsd[degs10,]),max.n.clus=15, n.clust=2, n.pca=4) #[degs10,]


#Use clus$grp to rename to in2in and off2off -
#Host
colnames(a.vsd) 
clus$grp=c(substr(colnames(a.vsd), start=1,stop=2)) #set the transplantation site as the groups


# now lets build a discriminant function for these two groups:
dp=dapc(t(a.vsd[degs10,]),clus$grp, n.da=1, perc.pca=80) #[degs10,]
# HOST: PCs: 6, functions: 1. For two groups only one discriminant function is possible.
# SYM: PCs: 6, functions: 1.

quartz()
scatter(dp,bg="white",scree.da=FALSE,legend=TRUE,solid=.4) #discriminant function for transplant


#APPLY FUNCTION TO OTHER SET OF SAMPLES
pred.sup<-predict.dapc(dp,newdata=(t(a.vsd.supp[degs10,]))) #skip IO11C for host b/c outlier sample in WGCNA

names(pred.sup)
pred.sup$assign
colnames(a.vsd.supp)

#must create another dataframe structure in order to plot these predicted values
test<-dp
test$ind.coord<-pred.sup$ind.scores
test$posterior<-pred.sup$posterior
test$assign<-pred.sup$assign


#HOST - for plotting distributions, must say which samples in which group
test$grp<-as.factor(substr(colnames(a.vsd.supp), start =1, stop=2)) #make sure origin is same num as before: IN=2, OFF=1


quartz()

scatter(test,bg="white",scree.da=FALSE,legend=TRUE,solid=.4,ylim=c(0,0.6),xlim=c(-4,4)) 
#adjust axes to correspond to previous graph, then overlay plots in photoshop





######## LOOK AT DISTRIBUTIONS ########
library(ggplot2)
t=pred.sup$ind.scores
h=dp$ind.coord
a=data.frame(rbind(t,h))
head(a)
a$treat=substr(rownames(a), start=1, stop=2)
mns=tapply(a$LD1, a$treat, mean)
#plot discriminant analysis
color.set=c('blue', 'turquoise', 'orange', 'red')
quartz()
ggplot(data=a, aes(LD1, fill=treat, color=treat)) + geom_density(alpha=0.6) + theme(panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +scale_color_manual(values=color.set) + scale_fill_manual(values=color.set)  #+ geom_vline(xintercept=mns, color=color.set)

######## SAVE RESULTS ########
ge.dapc = a

#for all genes
save(ge.dapc, file='~/gitreps/reciprocal_transplant_methylation/datasets/ge.dapc_all_genes.Rdata')

#for plastic genes
save(ge.dapc, file='~/gitreps/reciprocal_transplant_methylation/datasets/ge.dapc.Rdata')

#-------------------------------------------- Plot correlation with gain
lda.res=a
rownames(lda.res) = sub("_2m", "", rownames(lda.res))
#--------gather the mbd discim data ------------
lnames=load('bayRT_rlog_conditions_traits.RData')
lda.res$Colony.ID = rownames(lda.res)
head(lda.res)
x = merge(traits, lda.res, by='Colony.ID', all.x=T)
x=x[order(x$Colony.ID),]
sum(x$Colony.ID == traits$Colony.ID)==nrow(traits)
traits$ge.LD1=x$LD1
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

summary(lmkk)
summary(lmko)
summary(lmok)
summary(lmoo)
#########

head(traits)
ge.traits=traits
save(ge.traits, file='ge.traits.Rdata') 










#look at mortality
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
traits$oriMort=paste(paste(traits$ori, traits$tra, sep=""), traits$mortality, sep="_")
mtraits=traits[!is.na(traits$mortality),]
boxplot(mtraits$mbd.LD1~mtraits$oriMort)





head(traits)
traits$geno=paste(traits$ori, traits$num, sep='')
genos=unique(traits$geno)
g='K1'
g='O1'

otraits=traits[traits$ori=='O',]
for (g in genos){
	sub=otraits[otraits$geno==g,]
}

traits[traits$mortality==1,]


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













##retain DFA values for additional calculations
dpc=data.frame(rbind(dp$ind.coord,pred.sup$ind.scores))

write.csv(dpc,"Sym_IndividCoord_DFA_1174genes.csv",quote=F) #or Host_IndividCoord_DFA_7008genes.csv

###Testing significance of DFA differences - MCMCglmm

# coral 

coral=read.csv("Host_IndividCoord_DFA_7008genes.csv")

# setting up genotype (gt), origin (ori) and home/away (away) factors
coral$gt=gsub("i|o|A|C","",coral$X)
coral$ori=sub("oo|oi","o",coral$group)
coral$ori=sub("io|ii","i",coral$ori)
coral$away=0
coral$away[coral$group %in% c("oi","io")]=1

library(MCMCglmm)

# weak inverse wishart prior with parameter expansion for random effect of genotype (this is standard in MCMCglmm, the results are actually identical with the default uniform prior)
prior = list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V=1, nu=0.002,alpha.mu=0, alpha.V=1000)))
cd=MCMCglmm(LD1~ori+ori:away,random=~gt,data=coral,prior=prior,nitt=75000, thin=25, burnin=5000)
summary(cd)
            # post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# (Intercept)     4.395    3.738    5.057    942.0 <0.001 ***
# orio           -8.475   -9.388   -7.503    815.1 <0.001 ***
# orii:away      -6.115   -7.159   -5.150   1000.0 <0.001 ***
# orio:away       3.690    2.877    4.532    831.3 <0.001 ***

# calculating difference in magnitudes of orii:away and orio:away using sampled sets of parameters:
awayDelta=abs(cd$Sol[,"orii:away"])-abs(cd$Sol[,"orio:away"])

# 95% credible interval:
HPDinterval(awayDelta)
       # lower    upper
# var1 1.207062 3.717263

#MCMC p-value:
if (is.na(table(awayDelta<0)[2])) {
	cat("p <",signif(1/length(awayDelta),1))
} else { cat("p =",signif(table(awayDelta<0)[2]/length(awayDelta),2)) } # 3.6e-4

#------------------------------
# same for symbionts

sym=read.csv("Sym_IndividCoord_DFA_1174genes.csv")

sym$gt=gsub("i|o|A|C","",sym$X)
sym$ori=sub("oo|oi","o",sym$group)
sym$ori=sub("io|ii","i",sym$ori)
sym$away=0
sym$away[sym$group %in% c("oi","io")]=1

prior = list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V=1, nu=0.002,alpha.mu=0, alpha.V=1000)))
symd=MCMCglmm(LD1~ori+ori:away,random=~gt,prior=prior,data=sym,nitt=75000, thin=25, burnin=5000)
summary(symd)
            # post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# (Intercept)     2.575    2.012    3.168     1000 <0.001 ***
# orio           -4.699   -5.465   -3.866     1000 <0.001 ***
# orii:away      -2.756   -3.570   -1.963     1000 <0.001 ***
# orio:away       1.864    1.114    2.584     1000 <0.001 ***

# difference in magnitudes of orii:away and orio:away:
SymAwayDelta=abs(symd$Sol[,"orii:away"])-abs(symd$Sol[,"orio:away"])
HPDinterval(SymAwayDelta)
          # lower    upper
# var1 -0.2836332 1.953451

if (is.na(table(SymAwayDelta<0)[2])) {
	cat("p <",signif(1/length(SymAwayDelta),1))
} else { cat("p =",signif(table(SymAwayDelta <0)[2]/length(SymAwayDelta),2)) } # 0.069



############################################## OKAY, now ready for WGCNA 

datExpr0 = as.data.frame(t(dat[,1:45]))

gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK #if TRUE, no outlier genes

### Outlier detection incorporated into trait measures. 

traitData= read.csv("WGCNA_traitsHost.csv") 
dim(traitData)
head(traitData)
names(traitData)

#host
allTraits<-traitData


dim(allTraits)
names(allTraits)

# Form a data frame analogous to expression data that will hold the clinical traits.
Samples = rownames(datExpr0);
rownames(allTraits) = allTraits$Sample;
traitRows = match(Samples, allTraits$Sample);
datTraits = allTraits[traitRows, -1];


table(rownames(datTraits)==rownames(datExpr0)) #should return TRUE if datasets align correctly, otherwise your names are out of order


#sample dendrogram and trait heat map showing outliers
A=adjacency(t(datExpr0),type="signed")                 #SELECT SIGNED OR UNSIGNED HERE
# this calculates the whole network connectivity
k=as.numeric(apply(A,2,sum))-1
# standardized connectivity
Z.k=scale(k)
thresholdZ.k=-2.5 # often -2.5
outlierColor=ifelse(Z.k<thresholdZ.k,"red","black")
sampleTree = flashClust(as.dist(1-A), method = "average")
# Convert traits to a color representation where red indicates high values
traitColors=data.frame(numbers2colors(datTraits,signed=FALSE))
dimnames(traitColors)[[2]]=paste(names(datTraits))
datColors=data.frame(outlierC=outlierColor,traitColors)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree,groupLabels=names(datColors), colors=datColors,main="Sample dendrogram and trait heatmap")

# Remove outlying samples from expression and trait data - for host is 11C, zoox is 3A
remove.samples= Z.k<thresholdZ.k | is.na(Z.k)
datExpr=datExpr0[!remove.samples,]
datTraits=datTraits[!remove.samples,]

#file
save(datExpr, datTraits, file="symSamplesAndTraits1174signed_NOii3_Feb2016.RData") #OR hostSamplesAndTraits7008signed_NOoi11_Feb2016.RData

################Moving on!  Network construction and module detection


library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()
library(flashClust)
#disableWGCNAThreads()

#host
lnames = load(file="hostSamplesAndTraits7008signed_NOoi11_Feb2016.RData") 
#OR hostSamplesAndTraits7008signed_NOoi11_Feb2016.RData #symSamplesAndTraits1174signed_NOii3_Feb2016.RData

#Figure out proper SFT 

# # Choose a set of soft-thresholding powers
powers = c(seq(1,14,by=1)); #may need to adjust these power values to hone in on proper sft value
# # Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, networkType="signed", verbose = 2) #want smallest value, to plateau closest to 0.9 (but still under)

# # Plot the results:
windows()
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# # Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
    main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    labels=powers,cex=cex1,col="red");
# # this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# # Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
    xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
    main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

##########################

#####################IF too memory intensive for your CPU, run this code on supercomputer then skip to TOMstep below  #########
softPower=10 #For HOST: 6 ; for SYM: 10
adjacency=adjacency(datExpr, power=softPower,type="signed") #must change method type here too!!
#translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:
TOM= TOMsimilarity(adjacency,TOMType = "signed")
dissTOM= 1-TOM
######################open this file instead if adjacency and TOM run on supercomputer, otherwise skip to geneTree:
lnames = load(file="TOMstep_hostSamplesAndTraits7008signed_rlog_sft6_NOoi11_Feb2016.RData") #if you run the TOM and dissTOM steps elsewhere - load in files here
#######################################

geneTree= flashClust(as.dist(dissTOM), method="average")
#sizeGrWindow(10,6)
# pdf(file="dendrogram_thresh75_signed7k_11Cout_nov15.pdf", width=20, height=20)
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE,hang=0.04)
# dev.off()
#each leaf corresponds to a gene, branches grouping together densely are interconnected, highly co-expressed genes


minModuleSize=30 #we only want large modules
dynamicMods= cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize= minModuleSize)
table(dynamicMods)

#dynamicMods - HOST - bowtie2 - all 7k Genes based on 90% with less than 10% count rlog Feb16 - signed network - sft 6
# dynamicMods
   # 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 1781 1502 1311  814  252  209  195  188  178  155  116  106   63   63   44   31

#dynamicMods - SYM - bowtie2 - all 1100 Genes based on 90% with less than 10% count rlog Feb16 - signed network - sft 10

  # 1   2   3   4   5   6   7 
# 411 193 145 139 120  98  68 

dynamicColors= labels2colors(dynamicMods)
#plot dendrogram and colors underneath, pretty sweet
sizeGrWindow(8,6)
quartz()
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang= 0.05, main= "Gene dendrogram and module colors")

#Merge modules whose expression profiles are very similar or choose not to merge
#calculate eigengenes
MEList= moduleEigengenes(datExpr, colors= dynamicColors,softPower = 6) #change softPower for host or symbiont
MEs= MEList$eigengenes
#Calculate dissimilarity of module eigenegenes
MEDiss= 1-cor(MEs)
#Cluster module eigengenes
METree= flashClust(as.dist(MEDiss), method= "average")
#plot
quartz()
sizeGrWindow(7,6)
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")

MEDissThres= 0.15 #merge modules that are 85% similar
abline(h=MEDissThres, col="red")

merge= mergeCloseModules(datExpr, dynamicColors, cutHeight= MEDissThres, verbose =3)

mergedColors= merge$colors
mergedMEs= merge$newMEs

#pdf(file="HostMergeNetwork_sft8signed_Merge_11Cout_Nov15.pdf", width=20, height=20)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
#dev.off()

moduleColors= mergedColors
colorOrder= c("grey", standardColors(50))
moduleLabels= match(moduleColors, colorOrder)-1
MEs=mergedMEs

#save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file= "SymNetwork1k_rlog_signed_ii3out_Merge15_sft10_Feb16.RData")
#HostNetwork7k_rlog_signed_io11out_Merge15_sft6_Feb16.RData


#Relating modules to traits and finding important genes

#library(WGCNA)
# The following setting is important, do not omit.
#options(stringsAsFactors = FALSE);


# Load the expression and trait data saved in the first part
lnames = load(file = "hostSamplesAndTraits7008signed_NOoi11_Feb2016.RData"); #hostSamplesAndTraits7008signed_NOoi11_Feb2016.RData symSamplesAndTraits1174signed_NOii3_Feb2016.RData
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "HostNetwork7k_rlog_signed_io11out_Merge15_sft6_Feb16.RData"); #HostNetwork7k_rlog_signed_io11out_Merge15_sft6_Feb16.RData
#SymNetwork1k_rlog_signed_ii3out_Merge15_sft10_Feb16.RData
lnames


#######################Replot module dendrogram
quartz()
MEList= moduleEigengenes(datExpr, moduleColors, softPower=6)$eigengenes #change softPower for host=6 or symbiont=10
MEs<- MEList
#Calculate dissimilarity of module eigenegenes
MEDiss= 1-cor(MEs)
#Cluster module eigengenes
METree= flashClust(as.dist(MEDiss), method= "average")
#plot
sizeGrWindow(7,6)
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")

################################now for module trait heatmap
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

datTraits2=datTraits[,c(2:4,9:15)]

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors,softPower=6)$eigengenes #change softPower for host=6 or symbiont=10
MEs = orderMEs(MEs0)

#correlations of traits with eigengenes
moduleTraitCor = cor(MEs, datTraits2, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
Colors=sub("ME","",names(MEs))

#correlations of genes with eigengenes
moduleGeneCor=cor(MEs,datExpr)
moduleGenePvalue = corPvalueStudent(moduleGeneCor, nSamples);

quartz()
#represent module trait correlations as a heatmap
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                        signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
             xLabels = names(datTraits2),
             yLabels = names(MEs),
             ySymbols = names(MEs),
             colorLabels = FALSE,
             colors = blueWhiteRed(50),
             textMatrix = textMatrix,
             setStdMargins = FALSE,
             cex.text = 0.5,
             zlim = c(-1,1),
             main = paste("Module-trait relationships"))


#represent module trait correlations as a heatmap - version 2 = prettier 

# module-trait correlations
quartz()
library(RColorBrewer)
modLabels=sub("ME","",names(MEs))

ps=signif(moduleTraitPvalue,1)
cors=signif(moduleTraitCor,2)
textMatrix = cors;
#paste(cors, "\n(",ps, ")", sep = "");
textMatrix[ps>0.05]="-"
dim(textMatrix) = dim(moduleTraitCor)

par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
xLabels = names(datTraits2),
ySymbols = modLabels,
yLabels = modLabels,
colorLabels = FALSE,
colors = colorRampPalette(c("blue","lightblue","white","coral","red"))(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.7,
zlim = c(-0.7,0.7),
main = paste("Host Module-Trait correlations"))

# module size barplot
labelShift=300 # increase to move module size labels to the right
quartz()
par(mar = c(6, 8.5, 3, 3));
mct=table(moduleColors)
mct[modLabels]
x=barplot(mct[rev(modLabels)],horiz=T,las=1,xlim=c(0,4500),col=rev(modLabels))
text(mct[rev(modLabels)]+labelShift,y=x,mct[rev(modLabels)],cex=0.9) 

             
#Gene relationship to trait and important modules:

# Define variable weight containing the weight column of datTrait - leave weight as variable, but change names in first 2 commands
weight = as.data.frame(datTraits$Home1Away0); #change to your trait name of interest
names(weight) = "Home1Away0"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");

#Gene-trait significance correlation plots
#windows()

moduleCols=c( "brown","turquoise","blue") #HOST bowtie2 "midnightblue","greenyellow","black","turquoise"
#moduleCols =c("yellow",  "blue",  "brown", "grey", "turquoise"  ) #SYM bowtie2
par(mfrow=c(1,4))
for (module in moduleCols) {
column = match(module, modNames);
moduleGenes = moduleColors==module;
#sizeGrWindow(7, 7);
#par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
 abs(geneTraitSignificance[moduleGenes, 1]),
                 xlab = paste("ModMem in", module, "module"),
                 ylab = "Gene Sig for Home1Away0",
                 main = paste("MM vs. GS\n"),
                 cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
   }     

#plotting massive table of all information - module membership, genes, gene names, etc.

annot=read.table("kb8_mar2014_iso2gene.tab",header=FALSE,sep="\t") #past_apr2014_iso2gene.tab; kb8_mar2014_iso2gene.tab

probes = names(datExpr)
probes2annot = match(probes,annot$V1)
datGS.Traits=data.frame(cor(datExpr,datTraits,use="p"))
names(datGS.Traits)=paste("cor",names(datGS.Traits),sep=".")
datME=moduleEigengenes(datExpr,moduleColors)$eigengenes
datKME=signedKME(datExpr, datME, outputColumnName="MM.")
datOutput=data.frame(ProbeID=names(datExpr),annot[probes2annot,],moduleColors,datKME,datGS.Traits)
datOutput=datOutput[order((datOutput$MM.turquoise),decreasing=T),]
write.table(datOutput,"AnnotatedNetworkAnalysisResults1k_rlog_signed_ii3out_Merge15_sft10_Feb16.csv",row.names=F,sep=",") 


#####Alternatively, if coming back to file later
datOutput=read.csv("AnnotatedNetworkAnalysisResults1k_rlog_signed_ii3out_Merge15_sft10_Feb16.csv")
#AnnotatedNetworkAnalysisResults7k_rlog_signed_11Cout_Merge15_sft6_Feb16.csv
#AnnotatedNetworkAnalysisResults1k_rlog_signed_ii3out_Merge15_sft10_Feb16

#Creating files for GO analysis of modules
##############Execute this entire block of code through creating VSD files 
#Host: "pink","greenyellow","cyan","turquoise","red","black","yellow","midnightblue", "blue", "green", "lightcyan", "magenta", "purple", "salmon"
#Sym: "turquoise","blue","brown","yellow","red","black", "green"

Modcolors=c("pink","greenyellow","cyan") #change to vector of desired module color names
for (col in Modcolors) {
tab=datOutput[,c(1,4)]

tab2=tab   #rbind(tab,rest)

tab2$moduleColors=as.character(tab2$moduleColors)
#Categorical GeneOntology by Module
tab2$moduleColors[tab2$moduleColors!=col]<-0
tab2$moduleColors[tab2$moduleColors==col]<-1 
tab2$moduleColors=as.factor(tab2$moduleColors) 
print(col)
print(summary(tab2)) #do counts match table of module colors?
#tab2$moduleColors=as.numeric(tab2$moduleColors) 
print(head(tab2))

write.csv(tab2,file=paste("symGO_MM",col,"_categorical.csv", sep=""),quote=F,row.names=F) #change to host or sym

#Making VSD files by module for GO plot functions

vs=t(datExpr)
cands=names(datExpr[moduleColors==paste(col)])   

c.vsd=vs[rownames(vs) %in% cands,]
print(nrow(c.vsd)) #should correspond to module size
write.csv(c.vsd,file=paste("host_vsd_MM",col,".csv", sep=""),quote=F) #change to host or sym
}

############################################
###Visualization of Gene Networks


##############################heatmap of module expression with bar plot of eigengene, no resorting of samples...
#names(dis)
sizeGrWindow(8,7);
which.module="blue" #pick module of interest = those by transplant: pink, lightcyan, turquoise, yellow 
ME=MEs[, paste("ME",which.module, sep="")]
genes=datExpr[,moduleColors==which.module ] #replace where says subgene below to plot all rather than just subset

#quartz()
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(genes) ),nrgcols=30,rlabels=F, clabels=rownames(genes), rcols=which.module,)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="sample")


##########To output ME by sample

which.module="blue" #pick module of interest # black      blue     green   magenta      pink       red turquoise    yellow 
ME=MEs[, paste("ME",which.module, sep="")]

melightcyan<-ME #get ME for blue module above
meturq<-ME
meyellow<-ME
mepink<-ME
mecyan<-ME
meblack<-ME
mered<-ME
megreen<-ME
megreenyellow<-ME
mebrown<-ME
meblue<-ME

#host
meout<-data.frame(cbind(rownames(datExpr),megreen,meblack,megreenyellow,mered,mepink,meturq,meyellow,melightcyan,mecyan,meblue))

#Sym
meout<-data.frame(cbind(rownames(datExpr),meblack,mered,mebrown,meturq,meblue))

write.csv(meout,"MEbySample_symMods.csv",quote=F,row.names=F)

boxpls<-read.csv("MEbySample_hostMods.csv") #or symMods
boxpls$OriByTran=as.factor(boxpls$OriByTran)

boxplot(meblue~OriByTran,boxpls)

##############################heatmap of module expression with bar plot of trait of interest by sample...

sizeGrWindow(8,7);
which.module="greenyellow" #pick module of interest #black    blue   brown   green    magenta    pink     red
which.trait="ChlAperSA" #change trait of interest here
datTraits2=datTraits[order((datTraits$ChlAperSA),decreasing=F),]#change trait of interest here

trait=datTraits2[, paste(which.trait)]
genes=datExpr[,moduleColors==which.module ] #replace where says subgene below to plot all rather than just subset
genes=genes[rownames(datTraits2),]

#quartz()
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(genes) ),nrgcols=30,rlabels=F, clabels=rownames(genes), rcols=which.module,)
par(mar=c(5, 4.2, 0, 0.7))
barplot(trait, col=which.module, main="", cex.main=2,
        ylab="Chl A/cm2",xlab="sample")#change trait of interest here

################################ HOST correlation of ME by GO term vs Sym Density

MCs1=rep("turquoise",40) #cellular response to stress 40/88 top GO term
MCs6=rep("blue",17) #ribosomal subunit 17/23
MCs7=rep("blue",95) #for crs and ribo biogen combined

#Calculate eigengene for genes in "cellular response to stress" GO term (see go_plotresults script to pull these genes out of GO_MWU results) 
MEs0_CSR = moduleEigengenes(datExpr[,c("isogroup00216" , "isogroup00434" , "isogroup00534" , "isogroup00594" , "isogroup01439" , "isogroup01770" ,"isogroup02026",  "isogroup02116" , "isogroup02219" , "isogroup02270" , "isogroup02625" , "isogroup03196" , "isogroup03214" , "isogroup03438" , "isogroup03991" , "isogroup04098" , "isogroup05034" , "isogroup05085" , "isogroup07245" , "isogroup08305" , "isogroup08882" , "isogroup09304" , "isogroup09380" , "isogroup10226" , "isogroup10507" , "isogroup11784" , "isogroup12594" , "isogroup14814" , "isogroup18324" , "isogroup20143" , "isogroup20913" , "isogroup00099" , "isogroup00186" , "isogroup00530" , "isogroup02086" , "isogroup09887" , "GCKDGN102DNX2X" , "isogroup03592" , "isogroup08724" , "isogroup10024")], MCs1)$eigengenes

#Calculate eigengene for genes in "cellular response to stress" and "ribosome biogenesis" GO term (=ESR) 
MEs0_Comb = moduleEigengenes(datExpr[,c("isogroup00216" , "isogroup00434" , "isogroup00534" , "isogroup00594" , "isogroup01439" , "isogroup01770" ,"isogroup02026",  "isogroup02116" , "isogroup02219" , "isogroup02270" , "isogroup02625" , "isogroup03196" , "isogroup03214" , "isogroup03438" , "isogroup03991" , "isogroup04098" , "isogroup05034" , "isogroup05085" , "isogroup07245" , "isogroup08305" , "isogroup08882" , "isogroup09304" , "isogroup09380" , "isogroup10226" , "isogroup10507" , "isogroup11784" , "isogroup12594" , "isogroup14814" , "isogroup18324" , "isogroup20143" , "isogroup20913" , "isogroup00099" , "isogroup00186" , "isogroup00530" , "isogroup02086" , "isogroup09887" , "GCKDGN102DNX2X" , "isogroup03592" , "isogroup08724" , "isogroup10024","GE907583"  ,"isogroup00486", "isogroup02353", "isogroup04037" ,"isogroup04356", "isogroup05257", "isogroup05321","isogroup05419" ,"isogroup05585" ,"isogroup05747" ,"isogroup06151", "isogroup06875" ,"isogroup07626" ,"isogroup07722","isogroup07814", "isogroup08086" ,"isogroup08261", "isogroup08677", "isogroup08704" ,"isogroup09277" ,"isogroup09400","isogroup09430" ,"isogroup10142", "isogroup10312", "isogroup10510" ,"isogroup10806" ,"isogroup11071", "isogroup11122","isogroup11478", "isogroup11696", "isogroup11894", "isogroup12695" ,"isogroup12730" ,"isogroup12866" ,"isogroup13357","isogroup13643", "isogroup13958", "isogroup14126", "isogroup14272" ,"isogroup14591" ,"isogroup16307", "isogroup16358","isogroup16525", "isogroup16963", "isogroup17745", "isogroup18152", "isogroup18829", "isogroup19176", "isogroup21158","isogroup21216", "isogroup21643", "isogroup21922", "isogroup22099", "isogroup26119", "isogroup25942")], MCs7)$eigengenes

#Calculate eigengene for genes in "Small Ribosomal Subunit"
MEs0_SRS = moduleEigengenes(datExpr[,c("isogroup02353" ,"isogroup04356" ,"isogroup05257" ,"isogroup05585" ,"isogroup05747" ,"isogroup06108" ,"isogroup07814" ,"isogroup08261" ,"isogroup10510" ,"isogroup11071" ,"isogroup11478" ,"isogroup11696","isogroup12695" ,"isogroup14272" ,"isogroup16525" ,"isogroup21158" ,"isogroup21922") ], MCs6)$eigengenes


SymCorGO=cbind(rownames(datTraits),datTraits$Cells.cm2,MEs0_CSR,MEs0_Comb,MEs0_SRS)
names(SymCorGO)=c("Sample","SymDens","ME_CellRespStress","ME_Combined","ME_SmallRibSub")
SymCorGO$group<-c(rep("ii",11),rep("io",8),rep("oi",13),rep("oo",12))

quartz()
require(gridExtra)

p1<-ggplot(data=SymCorGO,aes(x=ME_CellRespStress,y=SymDens,color=group)) +
	geom_smooth(method = "lm", se=TRUE, color="black", formula=y~x)+
	geom_point() + theme_bw()
p4<-ggplot(data=SymCorGO,aes(x=ME_SmallRibSub,y=SymDens,color=group)) +
	geom_smooth(method = "lm", se=TRUE, color="black", formula=y~x)+
	geom_point() + theme_bw()

grid.arrange(p1,p4,ncol=2)

cor.test(SymCorGO$ME_CellRespStress, (SymCorGO$SymDens))
cor.test(SymCorGO$ME_SmallRibSub, (SymCorGO$SymDens))


###################DFA ESR genes only, to calculate plasticity

ESR=read.csv("Host_annotESRgenes_CRSriboBiogen.csv")
names<-data.frame(do.call('rbind', strsplit(as.character(ESR$X),'.',fixed=TRUE)))

ESR<-cbind(ESR,names)

rownames(ESR)<-ESR$X2

#dat=ESR[56:95,2:45] #56:95 is just CSR genes, n=40
#dat=ESR[1:55,2:45] #56:95 is just Ribosome genes, n=55

dat=ESR[,2:45]
nrow(dat) #95 annotated ESR genes
degs10<-rownames(dat)

a.vsd<-dat[,grep("A",colnames(dat))] #Corals in their native reefs
a.vsd.supp<-dat[,grep("C",colnames(dat))] #Foreign transplants, no XX for symbionts; no IO11C for corals

######################now use appropriate dataset for analysis
pcp=prcomp(t(a.vsd[degs10,]), retx=TRUE, center=TRUE, scale.=TRUE) #scale.=TRUE
scores=pcp$x
screeplot(pcp,bstick=T) # only 1st PC is non-random when using diff exp genes only; higher PCs non-random when using all data...

# adegenet: finding clusters (even though we know what clusters we want) - choose 10 PCs and 2 groups
clus=find.clusters(t(a.vsd[degs10,]),max.n.clus=15) #[degs10,]


#Use clus$grp to rename to in2in and off2off -
#Host 
clus$grp=c(2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1) #tell the DF which groups you want to cluster; in this case in2in and off2off

# now lets build a discriminant function for these two groups:
dp=dapc(t(a.vsd[degs10,]),clus$grp) #[degs10,]
# HOST: PCs: 10, functions: 1. For two groups only one discriminant function is possible.

quartz()
scatter(dp,bg="white",scree.da=FALSE,legend=TRUE,solid=.4) #discriminant function for ORIGIN type expression

#Now, add in transplants and see where they fall along continuum

#HOST
pred.sup<-predict.dapc(dp,newdata=(t(a.vsd.supp[degs10,]))) 


names(pred.sup)
pred.sup$assign
names(a.vsd.supp)

test<-dp
test$ind.coord<-pred.sup$ind.scores
test$posterior<-pred.sup$posterior
test$assign<-pred.sup$assign
#HOST
test$grp<-as.factor(c(2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1)) #make sure origin is same num as before: IN=2, OFF=1

quartz()

scatter(test,bg="white",scree.da=FALSE,legend=TRUE,solid=.4) #adjust axes to correspond to previous graph

##retain DFA values for additional calculations
dpc=data.frame(rbind(dp$ind.coord,pred.sup$ind.scores))

write.csv(dpc,"Host_IndividCoord_DFA_ESRGenes2.csv",quote=F)

################Now correlate difference in DFA of ESR genes (i.e. plasticity) with sym density

df=read.csv("Host_IndividCoord_DFA_CRSgenes.csv")
head(df)
ii=subset(df,group=="ii")
io=subset(df,group=="io")
oo=subset(df,group=="oo")
oi=subset(df,group=="oi")


#BUT does greater plasticity in ESR equal higher sym dens? yes, 
p1<-ggplot(data=df,aes(x=DAPC_diffESR,y=SymDens_at_IN,color=group)) +
	geom_smooth(method = "lm", se=TRUE, color="black", formula=y~x)+
	geom_point() + theme_bw()

cor.test(df$SymDens_at_IN, df$DAPC_diffESR)

p2<-ggplot(data=df,aes(x=DAPC_ESR_at_.offshore,y=SymDens_at_IN,color=group)) +
	geom_smooth(method = "lm", se=TRUE, color="black", formula=y~x)+
	geom_point() + theme_bw()

#what about genotypes with higher baseline expression? no
cor.test(df$SymDens_at_IN, df$DAPC_ESR_at_.offshore)

require(gridExtra)
grid.arrange(p1+theme(legend.position="none"),p2+theme(), ncol=2)



##############Calculate costs of plasticity
#relate fitness in a given environment to trait values in that environment and to plasticity between environments
#multiple regression: one predictor is trait value, other is difference
#partial regression coefficient for difference term provides estimate of plasticity cost while controlling for mean trait value
#negative indicates selection against plasticity (cost)
#positive indicates selection for plasticity

#From Hendry: We quantified the evidence for costs of canalization and plasticity by compiling all studies that have used the selection gradient method proposed by van Tienderen (1991; DeWitt etÂ al., 1998; Scheiner & Berrigan, 1998) or a derivative of it. At minimum, this method involves a multiple regression of relative fitness measured within environment e (we) against trait values expressed within eÂ (Xe) and trait plasticities across all environments (plX): 

#NOTE: relative fitness calculated as % gain after one year for individual divided by individual with maximum % gain within each environment - think of it as what the absolute best possible growth could be at inshore and offshore reefs (though max is similar)

plas=read.csv("PlasticityCosts.csv")
head(plas)
env1=subset(plas,environment=="inshore")
env2=subset(plas,environment=="offshore")
#FOR WEIGHT GAIN ~ Fitness in native reef - separate into inshore and offshore datasets
#relative fitness expressed as growth of individual relative to global max growth possible
#traits (GE data) are standardized (value-mean/SD) by population, so they are comparable (i.e. offshore mean negative, inshore mean positive)

m1=(lm((relFitnessNativeGainLocal)~StdTraitValue+stdDiffValue,env1))
summary(m1)
confint(m1,'stdDiffValue')

m2=(lm((relFitnessNativeGainLocal)~StdTraitValue+stdDiffValue,env2))
summary(m2)
confint(m2,'stdDiffValue')

#For Inshore pops => coeff DiffValue positive but NS=> Est:0.1295     Pval: 0.619  95%CI: -0.8252666 1.084181
#For Offshore pops => coeff DiffValue negative but NS Est:-0.02321    Pval:0.701    95%CI: -0.1604971 0.1140861
#Create dataframe for plotting outcomes
#######################################

