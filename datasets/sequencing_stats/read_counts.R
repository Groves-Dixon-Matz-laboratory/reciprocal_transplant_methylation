#read_counts.R
setwd("~/gitreps/reciprocal_transplant_methylation/datasets/sequencing_stats")
library(plotrix)

#------- look at read counts ----------
rawcounts = read.table("rawReadCounts.tsv")
colnames(rawcounts) = c('file', 'rawcount')

filtcounts=read.table("trimmed_read_counts.txt", header = T)
colnames(filtcounts) = c('file', 'filtcount')

filtcounts = read.table("filteredReadCounts.tsv")
colnames(filtcounts) = c('file', 'filtcount')



dat = merge(rawcounts, filtcounts, by = 'file')
head(dat)
nrow(dat)


filtcounts=read.table("trimmed_read_counts.txt", header = T)
colnames(filtcounts) = c('file', 'filtcount')

#plot histograms for all counts
par(mfrow=c(2,1))
hist(dat$rawcount,
	main='Raw Read Counts',
	xlab = 'Read Counts')
hist(dat$filtcount,
	main='Filtered Read Counts',
	xlab = 'Read Counts')

library(plotrix)
summary(dat$rawcount)
sum(dat$rawcount)
sd(dat$rawcount)
std.error(dat$rawcount)
summary(dat$filtcount)
sum(dat$filtcount)
std.error(dat$filtcount)

plot(filtcount~rawcount, data = dat)


#------- mapping efficiencies ----------
dig = read.table("adig_zoox_combo_mapping_efficiencies.txt")
mil = read.table("amil_zoox_combo_mapping_efficiencies.txt")
summary(dig)
summary(mil)
std.error(dig$V1)
std.error(mil$V1)



#------- duplicates ----------
mdups = read.table("Adig_zoox_combo_metRemovalMetrics.tsv")
udups=read.table("Adig_zoox_combo_ubRemovalMetrics.tsv")
dups = rbind(mdups, udups)
dim(dups)
colnames(dups) = c('sample', 'dup')
summary(dups)
std.error(dups$dup)