#protein_table.R
setwd("~/git_Repositories/reciprocal_transplant_methylation/")
ptable = read.table("datasets/ProteinTable10529_263537.txt", header = F, sep = "\t", quote = "")#note this file can be downloaded from here: https://www.ncbi.nlm.nih.gov/genome/proteins/10529?genome_assembly_id=263537
head(ptable)
ptable = ptable[,c(2,7,8,10)]
colnames(ptable) = c('contig', 'locusName', 'genbank.prot', 'protein.name')
head(ptable)
ptable=ptable[!duplicated(ptable$locusName),]    #remove lines for multiple isogroups
save(ptable, file='datasets/ProteinTable.Rdata')