#adegenet_snps.R
#run discriminate analysis of principal components on SNP data
#Groves Dixon
#last updated 8-28-17
#updated to use SNPs called with A.digitifera reference and mpileup

setwd("~/gitreps/reciprocal_transplant_methylation/")
library(vcfR)
library(adegenet)


#upload the vcf
# gll=vcfR2genlight(read.vcfR("datasets/recip_snps_4-28-17_noSingletons.recode.vcf"))#old one generated with GATK and A.mil reference. 
# gll=vcfR2genlight(read.vcfR("datasets/recipMeth_final_mindp20_maxMiss95.recode.vcf"))#New SNPs generated using mpileup and A.dig reference 8/28/17
gll=vcfR2genlight(read.vcfR("datasets/recipMeth_final_mindp5_maxMiss8.recode.vcf"))#New SNPs generated using mpileup and A.dig reference 8/28/17
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
col=c("dodgerblue", "firebrick")
s.class(pca$scores,pop(gll),col=col,axesell=F,cstar=0,grid=F)


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
save(snp.dp, file='datasets/snp.dapc.Rdata')

