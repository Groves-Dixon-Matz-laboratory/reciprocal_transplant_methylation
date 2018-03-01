#SAVE YOUR DIRECTORY NAME AND SET WD
library(scales)
library(plotrix)
directory<-"~/gitreps/reciprocal_transplant_methylation"
source("scripts/reciprocal_methylation_project_functions.R")
setwd(directory)


lnames = load("datasets/GE_3mo.RData")
lnames
ge.vsd = vsd
colnames(ge.vsd) = colnames(counts)
head(ge.vsd)



lnames=load("datasets/rld_PROMOTERS__p-1000_200.Rdata")
lnames
gbm.vsd = rld.df
names = sub("_2m", "", colnames(vsd))
colnames(gbm.vsd) = names
head(gbm.vsd)



#standardize the two datasets
s.in.both = colnames(ge.vsd)[colnames(ge.vsd) %in% colnames(gbm.vsd)]
g.in.both = rownames(ge.vsd)[rownames(ge.vsd) %in% rownames(gbm.vsd)]
ge.vsd=ge.vsd[g.in.both,s.in.both]
gbm.vsd=gbm.vsd[g.in.both,s.in.both]
dim(ge.vsd)
dim(gbm.vsd)


#reshape the data
library(reshape)
gbm <- melt(gbm.vsd, id=c(colnames(gbm.vsd)))
ge <- melt(ge.vsd, id=c(colnames(ge.vsd)))
colnames(gbm) = c('gene', 'sample', 'gbm')
colnames(ge) = c('gene', 'sample', 'ge')
#check that the rows are equivalent
sum(paste(gbm$gene, gbm$sample, sep="_") == paste(ge$gene, ge$sample, sep="_")) == nrow(gbm)
#assemble into single dataframe
d=gbm
d$ge = ge$ge

plot.lm = function(df, x, y){
	plot(df[,y]~df[,x], col='grey')
	lm1=lm(df[,y]~df[,x])
	abline(lm1, col='red')
	print(summary(lm1))
	return(lm1)
}


add.windows=function(df, n.windows){
	q=quantile(df$ge, seq(0, 1, 1/n.windows))
	gbm.mns=c()
	ge.mns=c()
	ses=c()
	for (i in 1:(length(q)-1)){
		left=q[i]
		right=q[i+1]
		sub = df[df$ge >= left & df$ge <= right,]
		gbm.mns =append(gbm.mns, mean(sub$gbm))
		ge.mns = append(ge.mns, mean(sub$ge))
		ses=append(ses, std.error(sub$gbm))
	}
	return(data.frame(gbm.mns, ge.mns, ses))
}


plot(d$gbm~d$ge, col = alpha('black', 0.01), pch=1, cex=0.75, xlim=c(-1, 10), ylim=c(-1, 10), xlab="Normalized Counts (Tag-seq)", ylab="Promoter Normalized Counts (MBD-seq)")
lm1=lm(d$gbm~d$ge)
summary(lm1)
text(x=8, y=10, "R2=0.00")


n.windows=20
mn.gbm = tapply(d$gbm, INDEX=d$gene, mean)
mn.ge = tapply(d$ge, INDEX=d$gene, mean)
mnd = data.frame(mn.gbm, mn.ge)
mnd=mnd[order(mnd$mn.ge),]
q=quantile(mnd$mn.ge, seq(0, 1, 1/n.windows))
ge= mnd$mn.ge
qge= mnd$mn.ge
for (i in 1:n.windows){
	left=as.numeric(q[i])
	right=as.numeric(q[i+1])
	qge[ge >= left & ge <= right]<-as.character(names(q[i]))
}
mnd$q=as.character(qge)
unique(as.character(mnd$q))
length(unique(as.character(mnd$q)))
boxplot(mnd$mn.gbm~mnd$q, outline=F)

od=d[order(d$ge),]
ge=od$ge
qge=od$ge
n.windows=20
q=quantile(ge, seq(0, 1, 1/n.windows))
q
for (i in 1:n.windows){
	left=q[i]
	right=q[i+1]
	qge[ge >= left & ge <= right]<-i
}
od$q=qge
head(od)
unique(od$q)
boxplot(od$gbm~od$q, notch=F, outline=F)

#----------- build an RPKM plot -----------#
#load the raw counts for the MBD-seq mapping
lnames=load('datasets/capturedVflowthroughResults_p-1000_200.Rdata')
mdat = data.frame(res)
head(ge.vsd)
mn.ge=data.frame(apply(ge.vsd, 1, mean))
colnames(mn.ge) = c('ge')
head(mdat)
m=merge(mdat, mn.ge, by = 0)
head(m)
plot(m$ge~m$log2FoldChange, col = alpha('black', 0.2))
lm1=lm(m$ge~m$log2FoldChange)
summary(lm1)
abline(lm1, col='red')

#------------ build the plots ----------------#
#all data points
plot(d$gbm~d$ge, col = alpha('black', 0.01), pch=1, cex=0.1, xlim=c(-1, 10), ylim=c(-1, 10), xlab="Normalized Counts (Tag-seq)", ylab="Normalized Counts (MBD-seq)", mgp=c(2.1,1,0))
lm1=lm(d$gbm~d$ge)
summary(lm1)
text(x=7.5, y=9.75, expression(paste("R"^"2", "=0.03****")))
#boxplots
boxplot(od$gbm~od$q, notch=F, outline=F, mgp=c(2.1,1,0), xlab="Transcription Quantiles", ylab="MBD-seq normalized Counts", axes=F)
axis(2, mgp=c(2.1,1,0));axis(1, at=1:20, labels=1:20, mgp=c(2.1,1,0));box()





#------------- build failure of within gene plots ----------------#
par(mfrow=c(1,2))
plot(d$gbm~d$ge, col = alpha('grey', 0.01), pch=1, cex=0.75, xlim=c(-1, 10), ylim=c(-1, 10), xlab="Normalized Counts (Tag-seq)", ylab="Normalized Counts (MBD-seq)", mgp=c(2.1,1,0))
abline(lm(d$gbm~d$ge), lwd=4, col='red')

genes = unique(as.character(d$gene))
length(genes)

corrs = c()
for (g in genes){
	sub=d[d$gene==g,]
	lms = lm(sub$gbm~sub$ge)
	abline(lms, lwd=0.3, col=alpha('black', 0.05))
	corrs = append(corrs, cor(sub$gbm, sub$ge))
}
abline(lm(d$gbm~d$ge), lwd=4, col='red')
mean(corrs)
boxplot(corrs)
abline(h=0, lty=2, col='grey')
text(x=1.3, y=.8, bquote("N="~.(length(corrs))))



#plastic gbm
pdat = d[d$gene %in% sig.genes,]
plot(d$gbm~d$ge)
lm1=lm(d$gbm~d$ge)
abline(lm1, col='red')
summary(lm1)

plot.lm(pdat, 'ge', 'gbm')
ws=add.windows(pdat, 20)
plotCI(ws$ge.mns, ws$gbm.mns, uiw=ws$ses, add=T, pch=26)


#pick gene set: 
bsig = sig.genes[sig.genes %in% pdat$gene]
genes = bsig
genes = unique(as.character(d$gene))
length(genes)

corrs = c()
for (g in genes){
	sub=d[d$gene==g,]
	lms = lm(sub$gbm~sub$ge)
	corrs = append(corrs, cor(sub$gbm, sub$ge))
}
mean(corrs)
boxplot(corrs)
abline(h=0, lty=2, col='grey')
text(x=1.4, y=.8, bquote("N="~.(length(corrs))))
corr.bak=corrs
save(corrs, file="~/Desktop/corrs.Rdata")

mean(corr.bak)
mn.gbm = tapply(pdat$gbm, pdat$sample, median)
mn.ge = tapply(pdat$ge, pdat$sample, median)
names(mn.gbm) == names(mn.ge)
# plot(mn.gbm~mn.ge)



t.gbm=gbm.vsd[rownames(gbm.vsd) %in% sig.genes,]
t.ge=ge.vsd[rownames(ge.vsd) %in% sig.genes,]

pc.gbm = prcomp(t.gbm)
pc.ge = prcomp(t.ge)
sgbm=pc.gbm$rotation
sge=pc.ge$rotation
rownames(sgbm) == rownames(sge)
treatment=substr(rownames(sgbm), start=1, stop=2)
cols = get.cols(treatment)
head(sgbm)
plot(sgbm[,1]~sge[,1], pch=21, bg=cols)





head(d)
num = substr(d$sample, start=3, stop=5)
ori = substr(d$sample, start=1, stop=1)
d$geno = paste(ori, num, sep="")
d$geneGeno=paste(d$gene, d$geno, sep="_")
d$treat=substr(d$sample, start=1, stop=2)
kk=d[d$treat=="KK",]
ko=d[d$treat=="KO",]
ok=d[d$treat=="OK",]
oo=d[d$treat=="OO",]

k = merge(kk, ko, by = 'geneGeno')
o = merge(oo, ok, by = 'geneGeno')
a = rbind(k, o)
a$gbmChange = a$gbm.x - a$gbm.y
a$geChange = a$ge.x - a$ge.y
abak=a

a=a[a$gene.x %in% sig.genes,]


plot(a$gbmChange~a$geChange)
lm1=lm(a$gbmChange~a$geChange)
summary(lm1)
abline(lm1, col='red')


mn.gbmc = tapply(a$gbmChange, INDEX=as.character(a$gene.x), mean)
length(mn.gbmc)
mn.gec = tapply(a$geChange, INDEX=as.character(a$gene.x), mean)

names(mn.gbmc) == names(mn.gec)

plot(mn.gbmc~mn.gec)
lm1=lm(mn.gbmc~mn.gec)
summary(lm1)
abline(lm1, col='red')
