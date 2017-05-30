#---------------- Upload the data ------------------
#SET UP THE DATA TO RUN DESEQ
library('DESeq2')
#SAVE YOUR DIRECTORY NAME AND SET WD
directory<-"~/git_Repositories/reciprocal_transplant_methylation/"
setwd(directory)

#READ IN THE COUNTS DATA
counts=read.table('datasets/all_rnaseq_counts.txt',header=TRUE,row.names=1) 
head(counts)
dim(counts)


#LOOK AT THE MEANS OF THE COLUMNS
#HOW MANY ARE GREATER THAN 3?
mns = apply(counts, 1, mean)
table(mns > 3)
# counts = counts[mns > 1,]
dim(counts)

#---------------- set up dataframe for sample varibles -----------------------
#BUILD A DATAFRAME ASSOCIATING SAMPLE NAMESWITH TREATMENT CONDITIONS
#set up sample names by clearing away .counts
sample = c()
for (i in colnames(counts)){
	sample = append(sample, strsplit(i, '.', fixed = T)[[1]][1])
}
sample

colony.id = paste(substr(sample, start=1, stop=1), substr(sample, start=3, stop = 10), sep = "")
data.frame(colnames(counts), colony.id)

#set up tranplant variable
transplant = sample
transplant[grep('KK', transplant)] <- 'K'
transplant[grep('OK', transplant)] <- 'K'
transplant[grep('KO', transplant)] <- 'O'
transplant[grep('OO', transplant)] <- 'O'
transplant

#set up origin variable
origin = sample
origin[grep('KK', origin)] <- 'K'
origin[grep('OK', origin)] <- 'O'
origin[grep('KO', origin)] <- 'K'
origin[grep('OO', origin)] <- 'O'
origin


treat=substr(sample, start=1, stop=2)

#build the dataframe with all the varialbes
COLDATA <- data.frame(sample, colony.id,transplant,origin, treat)
COLDATA

#####################

#change name to indicate this is captured only
countsco = counts

#BUILD A DESeq INPUT TABLE FROM THE DATAFRAME
ddsco<-DESeqDataSetFromMatrix(countsco,
	colData = COLDATA, 
	design = formula(~ transplant+origin))


# LRT testing
ddsco.o=DESeq(ddsco,test="LRT",reduced=~transplant)
ddsco.t=DESeq(ddsco,test="LRT",reduced=~origin)
orico=results(ddsco.o, contrast = c('origin', 'O', 'K'))
traco=results(ddsco.t, contrast=c('transplant', 'O','K'))
summary(orico)
summary(traco)
head(orico)
head(traco)

save(orico, file="datasets/orico_RNAseq.Rdata")
save(traco, file="datasets/traco_RNAseq.Rdata")



#####################################


### SET UP SUBSETS FOR DIFFERENTIAL EXPRESSION WITHIN TRANSPLANT SITES
#BUILD A DESeq INPUT TABLE FROM THE DATAFRAME
ddsco<-DESeqDataSetFromMatrix(countsco,
	colData = COLDATA, 
	design = formula(~ treat))

ddsco.o=DESeq(ddsco)
oriAtO=results(ddsco.o, contrast = c('treat', 'OO', 'KO')) #difference by origin at Orpheus
oriAtK=results(ddsco.o, contrast=c('treat', 'OK','KK'))    #difference by origin at Keppels
save(oriAtO, file="datasets/oriAtO_RNAseq.Rdata")
save(oriAtK, file="datasets/oriAtK_RNAseq.Rdata")

