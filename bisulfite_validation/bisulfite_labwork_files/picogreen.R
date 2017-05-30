# edit the  first line to match your input file name (2-column csv: well,fluorescence.)
# then mark and execute the whole script
# the output is printed out and saved under filename like input with prefix "res_"

#upload the picogreen fluorescence results
setwd("/Users/grovesdixon/git_Repositories/reciprocal_transplant_methylation/bisulfite_validation/bisulfite_labwork_files")
inFile="bs_amplicon_library_picogreen_correct.csv"
pg=read.csv(inFile, header=T)
names(pg)=c("well","fluo")
head(pg)

#assign the standard wells (B1-H1)
pg[,2]=pg[,2]-pg[pg$well=="A1",2]
pg=pg[pg$well!="A1",]
pg=pg[pg$fluo>0,]
pg$fluo=log(pg$fluo,10)
cal=pg[grep("[ABCDEFGH]1$",pg$well),]
cal$conc=log(c(100,100/3,100/3^2,100/3^3,100/3^4,100/3^5,100/3^6),10)
str(cal)

#plot the standard curve
plot(conc~fluo,cal,xaxt="n",yaxt="n",bty="n",xlab="fluorescence",ylab="DNA concentration (ng/ul)",mgp=c(2.1,1,0))
ll=lm(conc~fluo,cal)
abline(ll,col="red")
axis(1,at=cal$fluo,labels=round(10^(cal$fluo),0))
axis(2,at=cal$conc,labels=signif(10^(cal$conc),2))

#calculate sample concentrations from the curve
sam=pg[-grep("[ABCDEFGH]1$",pg$well),]
sam=sam[sam$fluo>0,]
sam$conc=predict(ll,newdata=data.frame(fluo=sam$fluo))
points(conc~fluo,sam,col="cyan3")
legend("topleft",col=c("black","cyan3"),pch=1,c("calibrarion","samples"),bty="n")
sam$conc=round(10^sam$conc,1)
sam$fluo=NULL
hist(sam[,2], xlab='ng/ul')
plot(density(sam[,2]),  xlab='ng/ul', main = '')

#output
names(sam)[2]="ng.ul"
sam
sam$row = substr(sam$well, start=1, stop=1)
sam$col = substr(sam$well, start=2, stop=6)
write.csv(sam,file=paste("res_",inFile,sep=""),quote=F,row.names=F)



#test correlation with previous messup run
sam.pg = pg[-grep("[ABCDEFGH]1$",pg$well),]

#upload other fluoresences
pg.bad = read.csv("./100ug_ml_standard_mistake/bs_amplicon_library_picogreen_mistake.csv", header = T)
names(pg.bad)=c("well","fluo")
head(pg.bad)
#assign the standard wells (B1-H1)
pg.bad[,2]=pg.bad[,2]-pg.bad[pg.bad$well=="A1",2]
pg.bad = pg.bad[pg.bad$well!="A1",]
pg.bad = pg.bad[pg.bad$fluo>0,]
pg.bad$fluo=log(pg.bad$fluo,10)
sam.pg.bad = pg.bad[-grep("[ABCDEFGH]1$", pg.bad$well),]


m.pg = merge(sam.pg, sam.pg.bad, by = 'well')
nrow(m.pg)
plot(m.pg$fluo.x~m.pg$fluo.y, pch=1)
lm1 = lm(m.pg$fluo.x~m.pg$fluo.y)
summary(lm1)
abline(lm1, col='red')



#check that concentrations match
resb = read.csv('./100ug_ml_standard_mistake/res_bs_amplicon_library_picogreen_mistake_standard_replaced.csv')
head(resb)
m.res = merge(sam, resb, by = 'well')
head(m.res)
nrow(m.res)
plot(m.res$ng.ul.x~m.res$ng.ul.y)
lm2=lm(m.res$ng.ul.x~m.res$ng.ul.y)
summary(lm2)
abline(lm2, col='red')

#So the results file res_bs_amplicon_library_picogreen_mistake_standard_replaced.csv is trustworthy




