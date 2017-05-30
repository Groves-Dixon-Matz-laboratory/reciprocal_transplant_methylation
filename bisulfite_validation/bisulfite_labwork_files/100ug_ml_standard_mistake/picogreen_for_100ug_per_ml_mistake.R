# edit the  first line to match your input file name (2-column csv: well,fluorescence. No header)
# then mark and execute the whole script
# the output is printed out and saved under filename like input with prefix "res_"
setwd("~/git_Repositories/reciprocal_transplant_methylation/bisulfite_validation/bisulfite_labwork_files/100ug_ml_standard_mistake")
inFile="bs_amplicon_library_picogreen_mistake.csv"

pg=read.csv(inFile,header=T)
names(pg)=c("well","fluo")
head(pg)
pg[,2]=pg[,2]-pg[pg$well=="A1",2]
pg=pg[pg$well!="A1",]
pg=pg[pg$fluo>0,]
pg$fluo=log(pg$fluo,10)
cal=pg[grep("[ABCDEFGH]1$",pg$well),]
low.keep = 3
cal=cal[(nrow(cal)-low.keep+1):nrow(cal), ]
cal
cal$conc=log(c(100,100/3,100/3^2,100/3^3,100/3^4,100/3^5,100/3^6),10)[(7-low.keep+1):7]
cal
# cal$conc=c(100,100/3,100/3^2,100/3^3,100/3^4,100/3^5,100/3^6)
str(cal)


ll=lm(conc~fluo,cal)
sam=pg[-grep("[ABCDEFGH]1$",pg$well),]
sam=sam[sam$fluo>0,]
sam$conc=predict(ll,newdata=data.frame(fluo=sam$fluo))

plot(conc~fluo,cal,xaxt="n",yaxt="n",bty="n",xlab="fluorescence",ylab="DNA concentration (ng/ul)",mgp=c(2.1,1,0), xlim=c(min(sam$fluo),max(cal$fluo)), ylim=c(min(sam$conc),max(cal$conc)), cex=2)
abline(ll,col="red")
xat=seq(from=min(sam$fluo), to=max(cal$fluo), length.out=6)
yat=seq(from=min(sam$conc), to=max(cal$conc), length.out=6)
axis(1,at=xat,labels=round(10^(xat),0))
axis(2,at= yat,round(10^yat,3)*50)
points(conc~fluo,sam,col="cyan3")
legend("topleft",col=c("black","cyan3"),pch=1,c("calibrarion","samples"),bty="n")


sam$conc=round(10^sam$conc,3)*50
sam$fluo=NULL
names(sam)[2]="ng/ul"
sam
write.csv(sam,file=paste("res_",inFile,sep=""),quote=F,row.names=F)
