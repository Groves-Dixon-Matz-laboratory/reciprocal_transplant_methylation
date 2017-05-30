#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================
# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = "~/git_Repositories/reciprocal_transplant_methylation";
setwd(workingDir); 
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);


# Load the expression and trait data saved in the first part
lnames = load(file = "wgcna/wgcna-01_output_p1000_200_iterate_mnCount10.RData");    
#The variable lnames contains the names of loaded variables.
lnames


# head(datTraits)
# dt2=read.csv("~/Desktop/change_zoox.csv")
# head(dt2)
# x=paste(dt2$genotype, '2m', sep="_")
# rownames(dt2) = x
# y=merge(datTraits, dt2, by = 0, all.x=T)
# dim(y)
# head(y)
# rownames(y) = y$Row.names
# sum(rownames(y)==rownames(datTraits))
# datTraits$zoox_change = y$change_zoox



# Load network data saved in the second part.
# lnames = load(file = "wgcna3b_manual_sft7_minModSize15_cutDpth0_signed.Rdata")
# lnames = load(file = "wgcna3b_manual_sft7_minModSize10_cutDpth0_signed.Rdata")
# lnames = load(file = "wgcna3b_manual_sft20_minModSize30_cutDpth0_signed.Rdata")
lnames = load(file = "wgcna/wgcna3b_manual_sft14_minModSize30_cutDpth0_signed.Rdata")
lnames = load(file = "wgcna/wgcna3b_manual_sft20_minModSize30_cutDpth0_signed.Rdata")
lnames

#match up the rownames in the expression 
rownames(datExpr)     #this tells us the rownames as given to WGCNA
rownames(datTraits)   #this is for the entire set of datTraits and includes samples that were not put into WGCNA
rownames(MEs) = rownames(datExpr)
datTraits = datTraits[rownames(datTraits) %in% rownames(datExpr),]  #reduce datTraits to those we have WGCNA results for
#if everything is matched up this should say TRUE
sum(rownames(datTraits) == rownames(datExpr)) == length(rownames(datExpr))
dim(datTraits)


#some final datTraits revising
datTraits$Colony.ID<-NULL
datTraits
time2zoox = datTraits$zoox
time2zoox[grep("3m", rownames(datTraits))]<-NA
time3zoox = datTraits$zoox
time3zoox[grep("2m", rownames(datTraits))]<-NA

datTraits$zoox<-NULL
datTraits$ZOOX=time2zoox
# datTraits$time3zoox=time3zoox

#make site-specific variables
#site-specific gain
gainInK = datTraits$gain
gainInO = datTraits$gain
gainInK[datTraits$transK == 0]<-NA
gainInO[datTraits$transK==1]<-NA
datTraits$gainInO = gainInO
datTraits$gainInK = gainInK
#site-specific mortality
datTraits$mortInK=datTraits$mortality
datTraits$mortInO=datTraits$mortality
datTraits$mortInK[datTraits$transK == 0]<-NA
datTraits$mortInO[datTraits$transK == 1]<-NA
datTraits$oriO = datTraits$oriK + 1
datTraits$oriO[datTraits$oriO == 2]<-0
datTraits = datTraits[,order(colnames(datTraits))]
datTraits


#add read counts
c = read.table("datasets/sequencing_stats/filteredReadCounts.tsv", header = F)
d = read.table("datasets/sequencing_stats/metRemovalMetrics_codingMapped.tsv")
head(d)
rownames(c) = c$V1
rownames(d) = d$V1
colnames(c) = c('sample', 'read_count')
colnames(d) = c('sample', 'dup_pct')
c$sample<-NULL
d$sample<-NULL
head(c)
head(d)
z = merge(c,d, by = 0, all = T)
rownames(z) = z$Row.names
z$Row.names<-NULL
head(z)
x = merge(z, datTraits, by = 0)
rownames(x) = x$Row.names
x$Row.names<-NULL
datTraits=x


remove = c("dup_pct", 'AMR', 'AMR.1', 'AMR.2', 'amr_t1.t2', 'amr_t1.t3', 'amr_t2.t3', 'AMR.1 AMR.2', 'area_t1', 'area_t2', 'area_t3')
datTraits = datTraits[,!colnames(datTraits) %in% remove]
datTraits



keep = c('gainInK', 'gainInO', 'mortInK', 'mortInO', 'oriK', 'oriO', 'ZOOX', 'zoox_change')
datTraits = datTraits[,colnames(datTraits) %in% keep]

#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate module eigengenes with color labels
#This will provide the first principal component for expression
#behavior for the genes in each module. On that principal component,
#each sample will have a loading value. These values can be corelated
#with our known sample traits to get an idea what biological mechanism
#the co-reguated genes in a given module might be responding to
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
#use the cor() function to get the correlations between the module eigengenes and the trait data
moduleTraitCor = cor(MEs, datTraits, use = "p");
#get p values as well
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


#write out module loadings
mEigs=MEs
rownames(mEigs) = rownames(datTraits)
# save(mEigs, file='datasets/moduleEigengenes.Rdata')

#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================

sizeGrWindow(10,6)

#now replot the heatmap
sizeGrWindow(10,6)

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                           signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(subCor)
par(mar = c(6, 8.5, 3, 3));

# Display the correlation values within a heatmap plot
COLORS = greenWhiteRed(50)
blueWhiteRed(50)


# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                           signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

# Display the correlation values within a heatmap plot
COLORS = greenWhiteRed(50)
blueWhiteRed(50)

#set up additional formatting variables
rows = rownames(moduleTraitCor)
sub.colors = substr(rownames(moduleTraitCor), 3, 50) #trim away the ME from the
module.sizes = paste(sub.colors, table(moduleColors)[sub.colors], sep = "\n")


labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = rownames(moduleTraitCor),
               ySymbols = module.sizes,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))



#--------plot the same thing only for ethanol significant modules--------
#first subset for the modules that are significant for ethanol 
CUT = 0.2
TRAIT = 'oriO'
TRAIT = 'transK'
TRAIT = 'ZOOX'
TRAIT = 'time2zoox'
TRAIT = 'gain'
TRAIT = 'gainInO'
TRAIT = 'gainInK'
TRAIT='mortInO'
TRAIT='mortality'
TRAIT='read_count'

#subset the correlations
subCor = moduleTraitCor[moduleTraitPvalue[, TRAIT] < CUT,]
subP = moduleTraitPvalue[moduleTraitPvalue[, TRAIT] < CUT,]
rows = rownames(moduleTraitCor)[moduleTraitPvalue[, TRAIT] < CUT]
sub.colors = substr(rownames(subCor), 3, 50) #trim away the ME from the
module.sizes = paste(sub.colors, table(moduleColors)[sub.colors], sep = "\n")

#now replot the heatmap
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(subCor, 2), "\n(",
                           signif(subP, 1), ")", sep = "");
dim(textMatrix) = dim(subCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
COLORS = greenWhiteRed(50)
blueWhiteRed(50)


labeledHeatmap(Matrix = subCor,
               xLabels = names(datTraits),
               yLabels = rownames(subCor),
               ySymbols = module.sizes,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
               
               
#=====================================================================================
#REPLOT THE MODULE EIGENGENE CLUSTERING TREE
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)

#plot them
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

#=====================================================================================

#=====================================================================================
#
#  Code chunk 4 Gather module membership data for the genes in each module
#
#=====================================================================================


# Define dataframe trait.df containing the a trait of interest from datTraits
trait.df = as.data.frame(datTraits[,TRAIT]);
names(trait.df) = TRAIT
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
head(geneModuleMembership)
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, trait.df, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(trait.df), sep="");
names(GSPvalue) = paste("p.GS.", names(trait.df), sep="");
modules = sub.colors
plot.cols = modules 

#=====================================================================================
#
#  Code chunk 5 plot scatterplots of module membership and trait correlations
#
#=====================================================================================   
#plot scatterplots for module
quartz()
par(mar=c(5,5))
m='salmon'
for (m in modules){
	column = match(m, modNames);
	moduleGenes = moduleColors==m;
	
	# sizeGrWindow(7, 7);
	# par(mfrow = c(1,1));
	verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
	                   abs(geneTraitSignificance[moduleGenes, 1]),
	                   xlab = paste("Module Membership in", m, "module"),
	                   ylab = paste("Correlation with", TRAIT),
	                   main = paste("Module membership vs. gene significance\n"),
	                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = 'black', bg = m, pch = 21, cex = 1.5)
}


#=====================================================================================
#
#  Code chunk 6 different outputting options
#
#=====================================================================================
#upload the protein annotations
lnames=load("datasets/ProteinTable.Rdata");lnames
head(ptable)






source("scripts/reciprocal_methylation_project_functions.R")
#get the results for a given module-trait combination
TRAIT='gainInO'
TRAIT='mortInO'
m='lightyellow'
tcors=get_module_corrs(m, TRAIT)
head(tcors, n=50)
head(tcors[order(tcors$trait_corr, decreasing=T),], n=50)
tcors[order(tcors$trait_corr, decreasing=T),] #look at top correlation genes






setwd("~/Desktop")
#output gene lists and build pdfs for normalized counts for module genes
genes = rownames(geneModuleMembership)
for (m in modules){
	moduleGenes = genes[moduleColors==m]
	print(paste("module =", m))
	print(length(moduleGenes))
	print(moduleGenes)
	outName = paste(m, 'genes.tsv', sep = "_")
	write.table(moduleGenes, outName, quote = F, row.names = F)
	# for (g in moduleGenes){
		# pdf(file = paste(paste(m, g, sep = "_"), "counts_figure.pdf", sep = "_"))
		# d <- plotCounts(dds.treatTime, gene=g, intgroup="treatTime", returnData=F)
		# dev.off()
	# }
}

#gene names
m = "turquoise"
genes = rownames(geneModuleMembership)
for (m in modules){
	moduleGenes = genes[moduleColors==m]
	write.table(moduleGenes, file = paste(m,'genes.txt', sep = "_"), quote = F, row.names = F)
}



#output for GO enrichment
lnames = load("all_annotations.Rdata")
head(all.anot)
all.anot = all.anot[!duplicated(all.anot[,'locusName']),]

# head(pdat)
# updat = pdat[!duplicated(pdat[,'locusName']),]
# head(updat)
# dim(pdat)
# dim(updat)

modules2 = c('red')
modules2 = c('lightyellow', 'magenta', 'yellow', 'midnightblue', 'blue', 'tan', 'green')
m = 'red'
for (m in modules2){
	print(paste("Module =", m))
	moduleGenes = moduleColors==m
	print(paste(sum(moduleGenes), "genes in module"))
	inMod = as.data.frame(geneModuleMembership[moduleGenes,])
	outMod = as.data.frame(geneModuleMembership[!moduleGenes,])
	inMem = data.frame(rownames(inMod), abs(inMod[,paste("MM", m, sep = "")]))
	print(paste("genes with membership data =", nrow(inMem)))
	colnames(inMem) = c('locusName', 'membership')
	outMem = data.frame(rownames(outMod), 0)
	colnames(outMem) = c('locusName', 'membership')
	res = rbind(inMem, outMem)
	res2 = merge(res, all.anot, by = 'locusName', all=F)[,c('prot.genbank', 'membership')]
	print(paste(nrow(res2[res2$membership >0,]), 'genes accounted for after merging'))
	head(res2)
	tail(res2)
	outName = paste(paste(m, TRAIT, sep ="_"), 'go_enrichment.csv', sep = "_")
	write.csv(res2, outName, quote = F, row.names = F)
}

head(geneModuleMembership)

outDir=""

m='salmon'
gnames = read.table("ProteinTable10529_263537.txt", header = F, sep = "\t")
gnames = gnames[!duplicated(gnames$locusName),]
colnames(gnames)[7]="locusName"
colnames(gnames)[10]="geneName"
head(gnames)
head()


#-------- PRINT OUT GENE NAMES ------------
m='turquoise'
m='salmon'
m='midnightblue'
moduleGenes = genes[moduleColors==m]
res=gnames[gnames$locusName %in% moduleGenes,]
unchar=res[grep("uncharacterized", res$geneName), 'geneName']
print(paste(length(unchar), "genes lack annotations"))
res2 = res[!res$geneName %in% (res[grep("uncharacterized", res$geneName), 'geneName']),]
print(paste(nrow(res2), 'with names:'))
res2

#=====================================================================================
#
#  Code chunk 6 plot export module data with gene names and descriptions
#
#=====================================================================================

#plot the counts for genes of interest in modules
#load the DESeq results
library('DESeq2')
lnames = load('/Users/grovesdixon/lab_files/coding4people/lafire/DESeq_all_counts_9-7-16/ethanol_LRT_results.Rdata') #load this for overall ethanol effect
lnames = load('/Users/grovesdixon/lab_files/coding4people/lafire/DESeq_all_counts_9-7-16/treatTime_LRT_results.Rdata')
lnames


#you can look at the actual difference in counts from the DESeq results for the module genes
modules
m = 'coral2'
moduleGenes = genes[moduleColors==m]
for (g in moduleGenes){
	# pdf(file = paste(g, "counts_figure.pdf", sep = "_"))
	d <- plotCounts(dds.eth, gene=g, intgroup="treatment", returnData=T)
	t = as.character(d$treatment)
	t[t == 'c'] <- 0
	t[t=='e'] <- 1
	d$treatment = as.numeric(t)
	plot(d$count~jitter(d$treatment, amount=.1), xlim = c(-.5,1.5))
	lm1 = lm(d$count~d$treatment)
	abline(lm1, col = 'red')
	
	# dev.off()
}
plot(d)



#for an individual gene of intersest
g = "ENSDARG00000055026"
plotCounts(dds.eth, gene=g, intgroup="treatment", returnData=F)
d <- plotCounts(dds.treatTime, gene='ENSDARG00000055026', intgroup="treatTime", returnData=F, main = 'ptch2 expression', axes = T, xlab = 'treatment groups')











#------------------- sanity check to see if ethanol modules are legitimate -------------------------
# permuts = 1000
# sanity_correlations = c()
# sanity_pvalues = c()

# for (i in 1:permuts){
	# sanity.eth = sample(datTraits$treatment)
	# datTraits$sanity.eth = sanity.eth
	# moduleTraitCor = cor(MEs, datTraits, use = "p");
	# moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
	# sanity_correlations = append(sanity_correlations, abs(moduleTraitCor[,'sanity.eth']))
	# sanity_pvalues = append(sanity_pvalues, moduleTraitPvalue[,'sanity.eth'])
# }

# #look at permuted correlations
# hist(sanity_correlations,
	# xlim = c(0, max(abs(subCor[,'treatment']))),
	# main = "Permuted Correlations")
# abline(v = abs(subCor[,'treatment']), col = 'red', lty = 2)

# #look at permuted p values
# hist(sanity_pvalues, main = 'Permuted P values')
# abline(v = abs(subP[,'treatment']), col = 'red', lty = 2)
# save(sanity_correlations, sanity_pvalues, file = 'permuted_module_stats.Rdata')


# #find permuted p values for the real ethanol module correlations
# subCor
# subP

# #function to iterate through set of correlations
# #and return the permutation value based on the sanity correlations for each
# getp = function(x){
	# res = c()
	# for (i in x){
		# y = length(sanity_correlations[sanity_correlations >= abs(i)]) / length(sanity_correlations)
		# print(paste('y=',y))
		# if (y == 0){
			# y = 1/length(sanity_correlations)
		# }
		# res = append(res, y)
	# }
	# return(res)
# }
# p = getp(subCor[,'treatment'])
# res = data.frame(rownames(subCor), p)
# res

#=====================================================================================























#======================================================================================
#
#
#
#======================================================================================
#output counts figures from DESeq analysis for the modules
#load DEseq results

lnames = load('/Users/grovesdixon/lab_files/coding4people/lafire/DESeq_all_counts_8-17-16/treatTime_LRT_results.Rdata')
lnames


#upload the gene ids
e2g = read.table("embl_to_gene.tsv", sep = "\t", header = T)
geneId = c()
for (i in e2g$embl_id){
	geneId = append(geneId, strsplit(i, ".", fixed = T)[[1]][1])
}
head(geneId)
e2g$embl = geneId
head(e2g)


#output gene lists and build pdfs for normalized counts for module genes
genes = rownames(geneModuleMembership)
for (m in modules){
	moduleGenes = genes[moduleColors==m]
	print(paste("module =", m))
	print(length(moduleGenes))
	outName = paste(m, 'genes.tsv', sep = "_")
	geneNames = e2g[e2g$embl %in% moduleGenes,]
	print('GeneNames:')
	print(geneNames)
	write.table(moduleGenes, outName, quote = F, row.names = F)
	# for (g in moduleGenes){
		# pdf(file = paste(paste(m, g, sep = "_"), "counts_figure.pdf", sep = "_"))
		# d <- plotCounts(dds.treatTime, gene=g, intgroup="treatTime", returnData=F)
		# dev.off()
	# }
}


              
#=====================================================================================
#
#  Code chunk 6
#   output data for Fisher tests using GO_MWU
#=====================================================================================


#choose modules of interest
# modules = c('bisque4', 'lightsteelblue', 'plum3')
genes = rownames(geneModuleMembership)



e2g = read.table("embl_to_gene.tsv", sep = "\t", header = T)
geneId = c()
for (i in e2g$embl_id){
	geneId = append(geneId, strsplit(i, ".", fixed = T)[[1]][1])
}
head(geneId)
e2g$embl = geneId

#funciton to output files for a set of modules
#for running Fisher's test with GO_MWU.R
isogroups = rownames(geneModuleMembership)
out_fisher_go = function(modules){
	for (m in modules){
		moduleGenes = moduleColors==m;
		goOut = data.frame(isogroups, as.numeric(moduleGenes))
		colnames(goOut) = c('isogroup', m)
		head(goOut)
		write.csv(goOut, paste(m, 'GO_input.csv', sep = "_"), quote = F, row.names = F)
	}
}
out_fisher_go(modules)



out_go = function(modules, out.type = 'normal'){
	for (m in modules){
		moduleGenes = moduleColors==m;
		mem.bool = as.numeric(moduleGenes)
		if (out.type == 'Fisher'){
			goOut = data.frame(isogroups, mem.bool)
			colnames(goOut) = c('isogroup', m)
			head(goOut)
			write.csv(goOut, paste(m, 'GO_Fisher_input.csv', sep = "_"), quote = F, row.names = F)
		}
		else {
			Fish_kME = mem.bool * geneModuleMembership[,paste("MM", m, sep = "")]
			goOut = data.frame(isogroups, Fish_kME)
			write.csv(goOut, paste(m, 'GO_kME_input.csv', sep = "_"), quote = F, row.names = F)
		}
	}
}
out_go(modules) #use this to output with 0 for nonmembers and kME for members
out_go(modules, out.type = 'Fisher') #use this for straight binary output for membership






names(datExpr)


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


names(datExpr)[moduleColors=="brown"]


#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================


annot = read.csv(file = "GeneAnnotation.csv");
dim(annot)
names(annot)
probes = names(datExpr)
probes2annot = match(probes, annot$substanceBXH)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.


#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================


# Create the starting data frame
geneInfo0 = data.frame(substanceBXH = probes,
                      geneSymbol = annot$gene_symbol[probes2annot],
                      LocusLinkID = annot$LocusLinkID[probes2annot],
                      moduleColor = moduleColors,
                      geneTraitSignificance,
                      GSPvalue)
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, weight, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight));
geneInfo = geneInfo0[geneOrder, ]


#=====================================================================================
#
#  Code chunk 10
#
#=====================================================================================


write.csv(geneInfo, file = "geneInfo.csv")



















