#snpRelate_recip_snps.R
#Groves Dixon


#see tutorial http://corearray.sourceforge.net/tutorials/SNPRelate/#data-analysis
##################################
####### Install Packages: ########
##################################
# source("http://bioconductor.org/biocLite.R")
# biocLite("gdsfmt")
# biocLite("SNPRelate")
##################################

#load packages
setwd("~/git_Repositories/reciprocal_transplant_methylation")
library(gdsfmt)
library(SNPRelate)

#reformat a vcf to run snpRelate
vcf.fn <-"datasets/recip_snps_1-27-17_noSingletons.recode.vcf"  #this one all all samples
snpgdsVCF2GDS(vcf.fn, "datasets/recip_snps_1-27-17_noSingletons.recode.vcf.gds", method="biallelic.only")
(genofile <- snpgdsOpen("datasets/recip_snps_1-27-17_noSingletons.recode.vcf.gds"))


#build the initial PCA data table and plot
pca <- snpgdsPCA(genofile, snp.id=NULL, num.thread=2, eigen.cnt = 4)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
tab <- data.frame(sample.id = pca$sample.id,
    EV1 = pca$eigenvect[,1],    # the first eigenvector
    EV2 = pca$eigenvect[,2],    # the second eigenvector
    EV3 = pca$eigenvect[,3],
    EV4 = pca$eigenvect[,4],
    stringsAsFactors = FALSE)
head(tab)
dim(tab)

#set up color coding
samples = tab$sample.id
origin = samples
origin[grep("K", origin)] <- 'dodgerblue'
origin[grep("O", origin)] <- 'firebrick'

#plot
plot(tab$EV1, tab$EV2, xlab="eigenvector 1", ylab="eigenvector 2", pch = 21, cex = 2, bg = origin, col='black')



#--------------------------------------------------------------------------

#####################################################
########### plot a color/symbol coded one ###########
#####################################################
color.set = c('darkblue', 'red', 'darkolivegreen')
color.set2x = c(color.set[1], color.set[1], color.set[2], color.set[2], color.set[3], color.set[3])
tab$site = 'none'
tab$site.col = 'black'
tab[grep('MCO', tab$sample.id), 'site'] <- 'Offshore'
tab[grep('MCD', tab$sample.id), 'site'] <- 'Deep'
tab[grep('MCN', tab$sample.id), 'site'] <- 'Inshore'
tab[grep('MCO', tab$sample.id), 'site.col'] <- color.set[3]
tab[grep('MCD', tab$sample.id), 'site.col'] <- color.set[1]
tab[grep('MCN', tab$sample.id), 'site.col'] <- color.set[2]

tab$age = 'none'
tab$age.shape = 0
age.pchs = c(22, 24)
legend.pchs = c(15, 17)
age.pchs2x = c(age.pchs[1], age.pchs[2], age.pchs[1], age.pchs[2], age.pchs[1], age.pchs[2])
tab[grep('A', tab$sample.id), 'age'] <- 'Adult'
tab[grep('J', tab$sample.id), 'age'] <- 'Juvenile'
tab[grep('A', tab$sample.id), 'age.shape'] <- age.pchs[1]
tab[grep('J', tab$sample.id), 'age.shape'] <- age.pchs[2]


#plot PCA
xlab = paste(paste('PCA1', round(pc.percent[1], 2), sep = "  "), "%", sep = "")
ylab = paste(paste('PCA2', round(pc.percent[2], 2), sep = "  "), "%", sep = "")
plot(jitter(tab$EV1, factor = 3000), tab$EV2, xlab=xlab, ylab=ylab, pch = as.numeric(tab$age.shape), bg = tab$site.col, cex = 3)
legend('bottomleft', c('deep adult', 'deep juvenile', 'inshore adult', 'inshore juvenile', 'offshore adult', 'offshore juvenile'), col = color.set2x, pch = legend.pchs)
#--------------------------------------------------------------------


#Use vegan to build a plot with elipses
# install.packages('vegan')
# install.packages('rgl')
# install.packages('ape')
library(vegan)
library(rgl)
library(ape)

tab$treat = paste(tab$site, tab$age, sep = " ")
deep = tab[tab$site == 'Deep',]
inshore = tab[tab$site == 'Inshore',]
offshore = tab[tab$site == 'Offshore',]
da = tab[tab$treat == 'Deep Adult',]
dj = tab[tab$treat == 'Deep Juvenile',]

na = tab[tab$treat == 'Inshore Adult',]
nj = tab[tab$treat == 'Inshore Juvenile',]

oa = tab[tab$treat == 'Offshore Adult',]
oj = tab[tab$treat == 'Offshore Juvenile',]

DRAW = 'lines'


#plot site PCA
xlab = paste(paste('PCA1', round(pc.percent[1], 2), sep = "  "), "%", sep = "")
ylab = paste(paste('PCA2', round(pc.percent[2], 2), sep = "  "), "%", sep = "")
scores = data.frame(tab$EV1, tab$EV2)

plot(tab$EV1, tab$EV2, xlab=xlab, ylab=ylab, pch = as.numeric(tab$age.shape), bg = tab$site.col, cex = 2)
# ordiellipse(tab[,2:3], tab$site, label=T, draw = c('lines'), col = color.set)
legend('bottomleft', c('deep adult', 'deep juvenile', 'inshore adult', 'inshore juvenile', 'offshore adult', 'offshore juvenile'), col = color.set2x, pch = legend.pchs)
ordiellipse(deep[,2:3], deep$site, label=F, draw = 'lines', col = color.set[1])
ordiellipse(inshore[,2:3], inshore$site, label=F, draw = 'lines', col = color.set[2], lwd = 2)
ordiellipse(offshore[,2:3], offshore $site, label=F, draw = 'lines', col = color.set[3], lwd = 2)
#----------------------------------------------------------------------------------------
######################################
############# run adonis #############
######################################
#this doesn't seem to work with missing values
#one solution is to use the random forest roughFix function, which just replaces missing values with modal allele for that locus
#create the original adat input files with vcf_2_pca_input.py
library(randomForest)


#upload the dataset that you want to run adonis on
adat = read.table("mcav_pca_input.tsv", header = T)                   #use this one for all samples
adat = read.table("mcav_adults_5-10-16_PCA_input.tsv", header = T)    #use this one for adults only
adat = read.table("mcav_juveniles_5-10-16_PCAinput.tsv", header = T)  #use this one for juveniles only


#rough fix the dataframe
head(adat)
adat = na.roughfix(adat)
head(adat)


#now transpose the data
adat = t(adat)


#set up condition data
site = row.names(adat)
site[grep('MCO', site)] <- 'Offshore'
site[grep('MCD', site)] <- 'Deep'
site[grep('MCN', site)] <- 'Inshore'

age = row.names(adat)
age[grep('A', age)] <- 'Adult'
age[grep('J', age)] <- 'Juvenile'

treat = paste(site, age, sep = "_")

condition = data.frame(site, age, treat)
head(condition)

#test for significant clustering
adonis(adat ~ site, data = condition, method = 'manhattan')


#test for sig clustering between inshore and offshore
head(condition)
head(adat)
no.deep = adat[condition$site != 'Deep',]
no.deep.cond = condition[condition$site != 'Deep',]
no.deep.cond$sanity = sample(no.deep.cond$site, length(no.deep.cond$site))
adonis(no.deep ~ site, data = no.deep.cond, method = 'manhattan')
adonis(no.deep ~ sanity, data = no.deep.cond, method = 'manhattan')





# dd.veg=vegdist(adat, "manhattan")
# div.dd.veg=dd.veg/1000
# head(div.dd.veg)

# dd.pcoa=pcoa(div.dd.veg) 
# head(dd.pcoa)
# scores=dd.pcoa$vectors
# quartz()
# plot(scores[,1], scores[,2],col=condition$site, xlab="PCo1", ylab="PCo2", main="PCoA by Reef")
############################################################################
############################################################################



################################################################
################## LOOK AT HABITAT MISMATCHES ##################
################################################################
#use cluster analysis on the PCs to assign groups
#if you want to run it for real comment this line in and continue
require(mclust)
sub = tab[tab$site != "Deep",]
x = tab$EV2
hist(x)

#use Bayesian Information Criterion to select the optimal mixture model/number of components
bic = mclustBIC(x, G = c(1:5), modelNames = c("E"))
summary(bic)
plot(bic)
#fit the model with two components
mod2 = Mclust(x, G = 2, modelNames = c('E'))
par(mfrow = c(2,1))
# plot(mod2, col = c('green', 'red'))
mod2$class


#isolate mismatched deep
head(tab)
d = tab[tab$site =='Deep',]
md = d[d$EV1 > -0.05,]
md[md$age == 'Adult',]
md[md$age == 'Juvenile',]
length(md$age[md$age == 'Adult'])
length(md$age[md$age == 'Juvenile'])

#isolate mismatched from inshore/offshore
head(tab)
i = tab[tab$site =='Inshore',]
o = tab[tab$site =='Offshore',]
mi = i[i$EV2 < 0,]
mo = o[o$EV2 > 0,]
mi[mi$age == 'Adult',]
mi[mi$age == 'Juvenile',]
mia = length(mi$age[mi$age == 'Adult'])
mij = length(mi$age[mi$age == 'Juvenile'])

mo[mo$age == 'Adult',]
mo[mo$age == 'Juvenile',]
moa = length(mo$age[mo$age == 'Adult'])
moj = length(mo$age[mo$age == 'Juvenile'])


#plot the number of habitat mismatches
par(mfrow = c(1,1))
m = as.table(cbind(c(mia, mij), c(moa, moj)))
colnames(m) = c('inshore', 'offshore'); rownames(m) = c("Adult", "Juvenile")
barplot(m, beside = T, names = c('Inshore', 'Offshore'))
par(xpd = T)
legend(4, 11, c('Adult', "Juvenile"), fill = grey.colors(2))



CUT = -0.05



#-------------- do fisher tests ----------------------

#function to compare number of juveniles and adults "mismatched" between habitats
do.fisher.test = function(site, cut, dir){
	i = tab[tab[,'site'] == site,]
	t = nrow(i)
	ta = nrow(i[i$age == 'Adult',])
	tj = nrow(i[i$age == "Juvenile",])
	if (dir == 'less'){
		mi = i[i$EV2 < cut,]
		}
	else{
		mi = i[i$EV2 > cut,]
		}
	# mi[mi$age == 'Adult',]
	# mi[mi$age == 'Juvenile',]
	mia = length(mi$age[mi$age == 'Adult'])
	mij = length(mi$age[mi$age == 'Juvenile'])
	mm = c(mia, mij)
	fit = c((ta-mia), c(tj - mij))
	in.tab = as.table(rbind(mm, fit))
	colnames(in.tab) = c('Adult', 'Juvenile')
	in.tab
	print(fisher.test(in.tab, alternative = 'less'))
	return(in.tab)
}
#run the function to look at age discrepancy for inshore and offshore
inshore = do.fisher.test('Inshore', 0, 'less')
offshore = do.fisher.test("Offshore", 0, 'greater')


#now combine the two together to see if, across all ishore and offshore there are more juvies
t.mm = c(10, 18)
t.fit = c(29, 20)
t.tab = as.table(rbind(t.mm, t.fit)); colnames(t.tab) = c('Adult', 'Juvie'); t.tab
barplot(t.mm, beside = T, col = grey.colors(2)); 
fisher.test(t.tab, alternative = 'less')
legend(.5, 15, c("Juveniles", "Adults"), fill = grey.colors(2))
#--------------------------------------------------------------------------
mns = mod2$parameters$mean
abline(v = mns[1], col = 'green', lwd = 3)
abline(v = mns[2], col = 'red', lwd = 3)

#for offshore
site.num = 1
site = 'Offshore'


#for inshore
site.num = 2
site = 'Inshore'

i = tab[tab[,'site'] == site,]
i$dist = abs(i$EV2 - mns[site.num])
ia = i[i$age == 'Adult',]
ij = i[i$age == 'Juvenile',]
boxplot(ia$dist, ij$dist, names = c('Adult', 'Juvenile'))
t.test(ia$dist, ij$dist, 'less')



ka = ia
kj = ij
ka = rbind(ka, ia)
kj = rbind(kj, ij)

boxplot(ka$dist, kj$dist, names = c('Adult', 'Juvenile'))
t.test(ka$dist, kj$dist, 'less')




#plot pairs of PCs up to 4
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(tab[,2:5], col=tab$site, pch = as.numeric(tab$age))


#Use vegan
# install.packages('vegan')
# install.packages('rgl')
# install.packages('ape')
library(vegan)
library(rgl)
library(ape)


tab$treat = paste(tab$site, tab$age, sep = " ")
deep = tab[tab$site == 'Deep',]
inshore = tab[tab$site == 'Inshore',]
offshore = tab[tab$site == 'Offshore',]
da = tab[tab$treat == 'Deep Adult',]
dj = tab[tab$treat == 'Deep Juvenile',]

na = tab[tab$treat == 'Inshore Adult',]
nj = tab[tab$treat == 'Inshore Juvenile',]

oa = tab[tab$treat == 'Offshore Adult',]
oj = tab[tab$treat == 'Offshore Juvenile',]

DRAW = 'lines'


#plot site PCA
xlab = paste(paste('PCA1', round(pc.percent[1], 2), sep = "  "), "%", sep = "")
ylab = paste(paste('PCA2', round(pc.percent[2], 2), sep = "  "), "%", sep = "")
scores = data.frame(tab$EV1, tab$EV2)

plot(tab$EV1, tab$EV2, xlab=xlab, ylab=ylab, pch = as.numeric(tab$age.shape), col = tab$site.col, cex = 1)
# ordiellipse(tab[,2:3], tab$site, label=T, draw = c('lines'), col = color.set)
legend('topleft', c('Deep Adult', 'Deep Juvenile', 'Inshore Adult', 'Inshore Juvenile', 'Offshore adult', 'Offshore juvenile'), col = color.set2x, pch = age.pchs2x)
ordiellipse(deep[,2:3], deep$site, label=T, draw = 'lines', col = color.set[1])
ordiellipse(inshore[,2:3], inshore$site, label=T, draw = 'lines', col = color.set[2], lwd = 2)
ordiellipse(offshore[,2:3], offshore $site, label=T, draw = 'lines', col = color.set[3], lwd = 2)




############################################################################
#
ibs <- snpgdsIBS(genofile, num.thread=2)
head(ibs)
summary(ibs)
image(ibs$ibs, col=terrain.colors(16))
ibs$sample.id

# install.packages('pheatmap')
library(pheatmap)
pheatmap(x, show_rownames = T)
quartz()

x = ibs$ibs
rownames(x) = ibs$sample.id
class(x)
heatmap(x, labRow = rownames(x), labCol = rownames(x))

############################################################################
#random forest
data(iris)
head(iris)
class(iris$Sepal.Length)
set.seed(71)    iris.rf <- randomForest(Species ~ ., data=iris, importance=TRUE, proximity=TRUE)


# install.packages("randomForest")
library(randomForest)
sdat = read.table("mcav_pca_input.tsv", header = T)
sdat = as.numeric(sdat)
head(sdat)
sdat = na.roughfix(sdat)
sdat = as.data.frame(t(sdat))
head(sdat)

#set up condition data
site = row.names(sdat)
site[grep('MCO', site)] <- 'Offshore'
site[grep('MCD', site)] <- 'Deep'
site[grep('MCN', site)] <- 'Inshore'

age = row.names(sdat)
age[grep('A', age)] <- 'Adult'
age[grep('J', age)] <- 'Juvenile'

treat = paste(site, age)

condition = data.frame(site, age, treat)
head(condition)


sdat$site = as.factor(site)
head(sdat)
rdat = sdat[,1:50]
class(sdat)
class(rdat$V1)
rdat$site = as.factor(site)
head(rdat)
class(rdat$site)
s.rf <- randomForest(site ~ ., data=rdat, importance=TRUE, proximity=TRUE) #don't run this on full dataset, it will freeze
round(importance(s.rf), 2)
head(s.rf$err.rate)
tail(s.rf$err.rate)
nrow(s.rf$err.rate)
ncol(rdat)


rdat = sdat[,1:500]
rdat$site = as.factor(site)
oob.record = c()
oob = 1
iter = 0
while  (oob > 0.05){
	iter = iter + 1
	s.rf <- randomForest(site ~ ., data=rdat, importance=TRUE, proximity=TRUE)
	oob = s.rf$err.rate[nrow(s.rf$err.rate), 'OOB']
	oob.record = append(oob.record, oob)
	imp = as.data.frame(round(importance(s.rf), 2))
	imp$max = apply(imp[,1:3], 1, max)
	keep = imp[imp$max > 0,]
	rdat = subset(rdat, select = rownames(keep))
	print(paste(paste('Round', iter), '...', sep = ""))
	print(paste("OBB =", oob))
	print(paste(ncol(rdat), 'variants remaining'))
	rdat$site = as.factor(site)
	if (iter > 1e2){
		print("Haven't reaced OO")
		break
	}
}

######################################
############# run adonis #############
######################################
#this doesn't seem to work with missing values
#not sure what to do about that


adat = read.table("mcav_beagled_pca_input.tsv", header = T) #create this file with vcf_2_pca_input.py
adat = adat[,2:length(adat)]
adat[adat=="none"] <- as.numeric(NA)
head(adat)
adat = t(adat)
head(adat)

# adat[adat=="none"] <- NA
# adat = as.numeric(adat)

#set up condition data
site = row.names(adat)
site[grep('MCO', site)] <- 'Offshore'
site[grep('MCD', site)] <- 'Deep'
site[grep('MCN', site)] <- 'Inshore'

age = row.names(adat)
age[grep('A', age)] <- 'Adult'
age[grep('J', age)] <- 'Juvenile'

treat = paste(site, age, sep = "_")

condition = data.frame(site, age, treat)
head(condition)

#test for significant clustering
adonis(adat ~ site, data = condition, method = 'manhattan')


#test for sig clustering between inshore and offshore
head(condition)
head(adat)
no.deep = adat[condition$site != 'Deep',]
no.deep.cond = condition[condition$site != 'Deep',]
no.deep.cond$sanity = sample(no.deep.cond$site, length(no.deep.cond$site))
adonis(no.deep ~ site, data = no.deep.cond, method = 'manhattan')
adonis(no.deep ~ sanity, data = no.deep.cond, method = 'manhattan')





dd.veg=vegdist(adat, "manhattan")
div.dd.veg=dd.veg/1000
head(div.dd.veg)

dd.pcoa=pcoa(div.dd.veg) 
head(dd.pcoa)
scores=dd.pcoa$vectors
quartz()
plot(scores[,1], scores[,2],col=condition$site, xlab="PCo1", ylab="PCo2", main="PCoA by Reef")
############################################################################
############################################################################










#plot full treatment PCA
scores = data.frame(tab$EV1, tab$EV2)
plot(tab$EV1, tab$EV2, xlab="PCA1", ylab="PCA2", pch = as.numeric(tab$age.shape), col = tab$site.col, cex = 1)
ordiellipse(tab[,2:3], tab$treat, label=T, draw = DRAW, col = color.set2x)
legend(-.2, -.02, c('Deep Adult', 'Deep Juvenile', 'Inshore Adult', 'Inshore Juvenile', 'Offshore adult', 'Offshore juvenile'), col = color.set2x, pch = age.pchs2x)

ordiellipse(da[,2:3], da$treat, label=T, draw = DRAW, col = color.set[1], lwd = 2)
ordiellipse(dj[,2:3], dj$treat, label=T, draw = DRAW, col = color.set[1], lwd = 2)
ordiellipse(oa[,2:3], oa$treat, label=T, draw = DRAW, col = color.set[3], lwd = 2)
ordiellipse(oj[,2:3], oj$treat, label=T, draw = DRAW, col = color.set[3], lwd = 2)
ordiellipse(na[,2:3], na$treat, label=T, draw = DRAW, col = color.set[2], lwd = 2)
ordiellipse(nj[,2:3], nj$treat, label=T, draw = DRAW, col = color.set[2], lwd = 2)



