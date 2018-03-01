#adegenet_snps.R
#run discriminate analysis of principal components on SNP data
#Groves Dixon
#last updated 8-28-17
#updated to use SNPs called with A.digitifera reference and mpileup

setwd("~/gitreps/reciprocal_transplant_methylation/")
library(vcfR)
library(adegenet)


#upload the vcf
gll=vcfR2genlight(read.vcfR("datasets/recipMeth_final_mindp5_maxMiss8_maf.2_maxNonref.8_FINAL.recode.vcf")); type = "ALL"
gll=vcfR2genlight(read.vcfR("datasets/genic_SNPs_maf.2_maxNonRef.8_maxMiss.8_mindp5_FINAL.recode.vcf")); type = 'GENIC'
class(gll)



#look at genlight object
x=as.matrix(gll)
gi=as.genind(x)


#plot the SNP index
# plot(gll)  #this takes a long time


#assign populations
pop=substr(gll@ind.names,1,1)
pop(gll)=pop
popnum=as.numeric(factor(pop,levels=c("O","K")))


######## PREPARE DATA FOR DAPC ########
pca=glPca(gll,nf=2)
quartz()
pca$scores[,1]=(-1)*pca$scores[,1]
col = get.cols(c("KO", 'OK'))
s.class(pca$scores,pop(gll),col=col,axesell=F,cstar=0,grid=F)
snp.pca = data.frame(pca$scores)

#plot pca
treat = substr(rownames(snp.pca), start = 1, stop=1)
treat[treat=="K"]<-"KO"
treat[treat=="O"]<-"OK"
variance <- pca$eig*100/sum(pca$eig)
plot(snp.pca$PC2~snp.pca$PC1, pch=21, bg=get.cols(treat), cex=1.5)




# adegenet: finding clusters (even though we know what clusters we want) - choose 4 PCs and 2 groups
clus=find.clusters(gll,max.n.clus=15, n.clust=2, n.pca=4)
clus$grp=pop

######## BUILD DISCIMINATE FUNCTION ########
dp=dapc(gi, pop=pop, n.da=1, perc.pca=80) 
#discriminate between KK and OO
#keep 1 discriminate function
#use PCs to account for 80% of var

#plot distributions
quartz()
scatter(dp,bg="white",scree.da=FALSE,legend=TRUE,solid=.4)

#plot my way
#plot with vanilla
library(scales)
a=data.frame(dp$ind.coord)
a$pop=pop
ldens=tapply(a$LD1, a$pop, density)
allx <- unlist(lapply(ldens, function(e) e$x))
ally <- unlist(lapply(ldens, function(e) e$y))
grp=as.factor(a$pop)
levels(grp)
MGP=c(2.1, .75, 0.0)
plot(allx, ally, type = "n", xlab ="Discriminant function", ylab = "Density", mgp=MGP, axes=F)
axis(1, mgp=MGP);axis(2, mgp=MGP)
colors=a$pop
colors[colors=='K']<-'blue'
colors[colors=='O']<-'red'
color.set=c('blue', 'red')
for (i in 1:length(ldens)) {
	polygon(c(ldens[[i]]$x, rev(ldens[[i]]$x)), c(ldens[[i]]$y, 
	                  rep(0, length(ldens[[i]]$x))), col = alpha(color.set[i], 0.6), 
	                  lwd = 1, border = color.set[i])
}
mns = tapply(a$LD1, a$pop, mean)
abline(v=mns)
snp.dp=dp
out=paste(paste("datasets/snp.dapc", type, sep="."), "Rdata", sep=".")
save(snp.dp, snp.pca, file=out)

