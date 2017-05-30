#Plot Tranche Results

#set the working directory
setwd("~/git_Repositories/reciprocal_transplant_methylation")

#upload the tranch data
dat = read.csv("./datasets/gatk_outputs/recalibrate_SNP.tranches.csv", header = T)
head(dat)

#plot the barplot
barplot(rev(dat$numNovel),
 horiz = T,
 xlab = "Number of Novel Variants",
 names = paste(rev(dat$targetTruthSensitivity), "%", sep = ""),
 main = "Novel Variants Aquired with Each Tranch",
 ylab = "Percentage of Known Variants Included",
)

 
 
 
#plot the line graph
quartz()
plot(dat$novelTiTv~dat$targetTruthSensitivity,
	ylab = "Ti/Tv for Novel Variants",
	xlab = "Percentage of Known Variants Included")
lines(dat$novelTiTv~dat$targetTruthSensitivity)