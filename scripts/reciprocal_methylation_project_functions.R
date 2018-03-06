#reciprocal_methylation_project_functions.R

#some global plotting variables
MGP=c(2.1, .75, 0.0)
color.set=c('blue', 'deepskyblue', 'darkgoldenrod1', 'firebrick1')
CEX.AXIS=1.5
# MGP=c(2.1,1,0)    #misha's


volcano_plot = function(dat, pcol='pvalue', log2col='log2FoldChange', fdrcol='padj', cut=0.1, XLIM=NULL, YLIM=NULL, MAIN='', draw.box=T, LEGEND='FDR<0.1'){
	dat=data.frame(dat)
	plot(-log(dat[,pcol], 10)~dat[,log2col], cex=0.5, mgp=MGP, xlab=expression(paste("Log"[2], ' Fold Difference Methylation')), ylab=expression(paste('log'[10],'(p-value)')), axes=F, xlim=XLIM, ylim=YLIM, main=MAIN)
	axis(1, mgp=MGP);axis(2, )
	if(draw.box){box()}
	sub=dat[!is.na(dat[,fdrcol]),]
	sub2=sub[sub[,fdrcol]<cut,]
	points(-log(sub2[,pcol], 10)~sub2[,log2col], pch=21,col="black",cex=0.6, bg='red')
	legend('topleft', LEGEND, pt.bg='red', col='black', pch=21, pt.cex=1.2, box.lwd=0, inset=c(0, -.2), xpd=T)
}



my_boxplot = function(df, ycol, index, YLAB='Match'){
	df[,index] = as.factor(df[,index])
	# levels(df[,index]) = c("KK", "OK", "KO", "OO")
	b=boxplot(df[,ycol]~df[,index], ylab=YLAB, mgp=MGP, outline=T, pch='', cex.axis= CEX.AXIS, cex.lab=CEX.AXIS)
	colors =get.cols(as.character(df[,'treat']))
	uni.treats = as.factor(unique(df[,'treat']))
	xs=as.numeric(as.factor(df[,index]))
	points(df[,ycol]~jitter(xs), pch=21, bg=colors, cex=1.5)
}


emmix = function(x, means,sigma,lambda,k){
  mod <- normalmixEM(x, mu = means, sigma = sigma, lambda = lambda, k=k, arbvar=T)
  summary(mod)
  par(xpd = T)
  plot(mod, which = 2, breaks = 30, density = TRUE, title = "\n") 
  return(mod)
}


# lnames = load("/Users/grovesdixon/lab_files/projects/recip_meth/deseq_11-4-16/dig_combo_mapped/mbd_classes.Rdata")
lnames = load("/Users/grovesdixon/lab_files/projects/recip_meth/deseq_12_13_16/mbd_classes_p-1000_200.Rdata")
plot_subset_mbd_histogram=function(subset, color, BREAKS=50){
	x = hist(classes$mbd.score, breaks = 70, xlim = c(-7, 7),main = "", xlab = 'Log2 Fold Difference')
	sub=classes[classes$gene %in% subset,]
	hist(sub$mbd.score, col=color, add=T, breaks=BREAKS)
}

plot_subset_mbd_density=function(subset, color, YLIM=c(0, 0.2), MAIN = 'Density', XLAB='Absolute Methylation'){
	sub=classes[classes$gene %in% subset,]
	plot(density(sub$mbd.score), col=color, ylim = YLIM, main = MAIN, lwd=1, xlab = XLAB)
	lines(density(classes$mbd.score), lwd=2)
	lines(density(sub$mbd.score), col=color, lwd=2)
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}



#funciton to plot read counts and add gene names for a set of genes
#uses the protein annotation table above to get protein product names
named_plotcounts=function(ddsco, sig, INTGROUP, LINE=T, returnData=T){
	s=data.frame(sig)
	s$locusName=rownames(s)
	s=s[order(s$pvalue, decreasing=T),]
	z=merge(s,ptable,by="locusName")
	z=z[!duplicated(z$locusName),]
	dif=nrow(s)-nrow(z)
	print(paste(nrow(s), 'genes in subset'))
	print(paste(nrow(z), 'have annotations'))
	plotCount=0
	newPlots=c()
	if(row(s) > 9){newPlots=seq(from=10, to=nrow(s), by = 10)}
	print(paste('quartz calls to be made =', length(newPlots)+1))
	quartz()
	for (g in rownames(s)){
		plotCount=plotCount + 1
		p=plotCounts(ddsco, gene=g, intgroup=INTGROUP, returnData=T)
		plot(p$count~jitter(as.numeric(p[,INTGROUP]), amount=.05),xlim=c(.5,2.5), axes=F, xlab=INTGROUP,ylab='normalized counts',main=g)
		axis(1, at=c(1,2), labels=levels(p[,INTGROUP]))
		axis(2)
		lm1=lm(p$count~as.numeric(p[,INTGROUP]))
		if(LINE){abline(lm1, col='red')}
		geneName=ptable[ptable$locusName==g,'protein.name'][1]
		title(main=paste("\n\n", as.character(geneName)))
		pvalue=round(s[g,'padj'], digits=3);title(xlab=paste("\nadjp=", pvalue), adj=1)
		pvalue=round(s[g,'pvalue'], digits=3);title(xlab=paste("p=", pvalue), adj=1, line = 1)
	}
	if(returnData){
		return(p)
	}
}

named_plotcounts_multipanel=function(ddsco, sig, INTGROUP, LINE=T){
	s=data.frame(sig)
	s$locusName=rownames(s)
	s=s[order(s$pvalue, decreasing=T),]
	z=merge(s,ptable,by="locusName")
	z=z[!duplicated(z$locusName),]
	dif=nrow(s)-nrow(z)
	print(paste(nrow(s), 'genes in subset'))
	print(paste(nrow(z), 'have annotations'))
	plotCount=0
	newPlots=seq(from=1, to=nrow(s), by = 9)
	print(paste('quartz calls to be made =', length(newPlots)+1))
	print("Annotated Gene Set:")
	z=z[order(z$pvalue),]
	print(z)
	for (g in rownames(s)){
		plotCount=plotCount + 1
		if(plotCount %in% newPlots){quartz();par(mfrow=c(3,3))}
		p=plotCounts(ddsco, gene=g, intgroup=INTGROUP, returnData=T)
		plot(p$count~jitter(as.numeric(p[,INTGROUP]), amount=.05),xlim=c(.5,2.5), axes=F, xlab=INTGROUP,ylab='normalized counts',main=g)
		axis(1, at=c(1,2), labels=levels(p[,INTGROUP]))
		axis(2)
		lm1=lm(p$count~as.numeric(p[,INTGROUP]))
		if(LINE){abline(lm1, col='red')}
		geneName=ptable[ptable$locusName==g,'protein.name'][1]
		title(main=paste("\n\n", as.character(geneName)))
		pvalue=round(s[g,'pvalue'], digits=3);title(xlab=paste("p=", pvalue), adj=1, line = 1)
		pvalue=round(s[g,'padj'], digits=3);title(xlab=paste("\nadjp=", pvalue), adj=1)
	}
}







plot_meth_rna_scatterplots=function(odat, CUT, P.TYPE, DO.BOTH=T){
	#origin
	par(mfrow=c(1,4))
	sub2 = odat[odat[,paste('rna', P.TYPE, sep="_")] < CUT,]
	plot(sub2[,'rna_log2FoldChange'] ~ sub2[,'met_log2FoldChange'], pch = 8, cex = 0.8, col = 'red')
	lm1 = lm(sub2[,'rna_log2FoldChange'] ~ sub2[,'met_log2FoldChange'])
	abline(lm1, col = 'red');abline(h=0,v=0,lty=2, col='grey')
	r= cor.test(sub2[,'rna_log2FoldChange'], sub2[,'met_log2FoldChange'], method = 'spearman')
	sig.rna = c('sig.rna', r$estimate, r$p.value)
	
	#subset for significant meth
	sub3 = odat[odat[,paste('met', P.TYPE,sep="_")] < CUT,]
	plot(sub3[,'rna_log2FoldChange'] ~ sub3[,'met_log2FoldChange'], pch = 8, cex = 0.8, col = 'blue')
	lm1 = lm(sub3[,'rna_log2FoldChange'] ~ sub3[,'met_log2FoldChange'])
	abline(lm1, col = 'blue');abline(h=0,v=0,lty=2, col='grey')
	
	#plot for all datapoints
	do.plot(odat, 'rna_log2FoldChange', 'met_log2FoldChange', "log2 Expression Diff", "log2 Methylation Diff", ylim = c(-5,5), main = 'Origin', xlim=c(-2,2))
	r=cor.test(odat[,'rna_log2FoldChange'], odat[,'met_log2FoldChange'], method = 'spearman')
	all=as.character(c('all.data', r$estimate, r$p.value))
	points(sub2[,'rna_log2FoldChange'] ~ sub2[,'met_log2FoldChange'], pch = 8, cex = 0.8, col = 'red')
	lm1 = lm(sub2[,'rna_log2FoldChange'] ~ sub2[,'met_log2FoldChange'])
	abline(lm1, col = 'red')
	r=cor.test(sub3[,'rna_log2FoldChange'], sub3[,'met_log2FoldChange'], method = 'spearman')
	sig.meth=as.character(c('sig.meth', r$estimate, r$p.value))
	points(sub3[,'rna_log2FoldChange'] ~ sub3[,'met_log2FoldChange'], pch = 8, cex = 0.8, col = 'blue')
	lm1 = lm(sub3[,'rna_log2FoldChange'] ~ sub3[,'met_log2FoldChange'])
	abline(lm1, col = 'black')
	
	results=rbind(all, sig.rna, sig.meth)
	colnames(results) = c('subset', 'rho', 'p.value')
	
	#subset for significant in both
	if (DO.BOTH == TRUE){
		sub4 = sub3[sub3[,paste('rna', P.TYPE, sep="_")] < CUT,]
		points(sub4[,'rna_log2FoldChange'] ~ sub4[,'met_log2FoldChange'], pch = 19, cex = 0.8, col = 'purple')
		lm1 = lm(sub4[,'rna_log2FoldChange'] ~ sub4[,'met_log2FoldChange'])
		abline(lm1, col = 'purple')
		r.both=cor.test(sub4[,'rna_log2FoldChange'], sub4[,'met_log2FoldChange'], method = 'spearman')
		sig.both=c('sig.both', r.both$estimate, r.both$p.value)
		do.plot(sub4, 'rna_log2FoldChange', 'met_log2FoldChange', "log2 Expression Diff", "log2 Methylation Diff", ylim = c(-5,5), main = 'Origin');abline(h=0,v=0,lty=2, col='grey')
		points(sub4[,'rna_log2FoldChange']~sub4[,'met_log2FoldChange'], pch=19, cex=1.2, col='purple')
		
		#assemble stats for each 
		results = rbind(all, sig.rna, sig.meth, sig.both)
		colnames(results) = c('subset', 'rho', 'p.value')
		return(results)
	}
	else{
		results=rbind(all, sig.rna, sig.meth)
		colnames(results) = c('subset', 'rho', 'p.value')
		return(results)
	}
}


remove_min_lines=function(dat, minimum){
	print("Dimentions before:")
	print(dim(dat))
	z = apply(dat, 1, sum)
	res = dat[z > minimum,]
	print("Dimentions after removal:")
	print(dim(res))
	return(res)
}



#plot correlations
do.plot = function(dat, y, x, ylab, xlab, xlim=NULL, ylim = NULL, main = '\n'){
	plot(dat[,y]~dat[,x],
	xlab = xlab,
	ylab = ylab,
	ylim = ylim,
	xlim = xlim,
	main = main,
	col = 'grey50')
	lm1 = lm(dat[,y]~dat[,x])
	print(summary(lm1))
	print(cor.test(dat[,y], dat[,x], method = 'spearman'))
	abline(lm1, col = 'grey')
}




#function to test if counts of significance are correlated between two rows of a dataframe
#Arguments:
#dataIn = the dataframe
#col1 = string for column 1 to look for significance in
#col2 = string for column 2 to look for significance in
#CUT = significance cutoff
#shared.col = a column in the dataframe with identifiers (this can actually be any column)
sig_overlap = function(dataIn, col1, col2, CUT, shared.col){
	dat= dataIn[,c(shared.col, col1, col2)]
	#get total
	all=dat[!is.na(dat[,col1]),]
	all=all[!is.na(all[,col2]),]
	
	#isolate significant rows for col1
	sig1 = dat[!is.na(dat[,col1]), ]
	sig1 = sig1[sig1[,col1] <= CUT,]
	
	#isolate non-significant rows for col1 (include NAs)
	notsig1 = dat[dat[,col1] > CUT,]
	
	
	#check we have all genes
	sum(nrow(sig1), nrow(notsig1))
	nrow(dat)
	
	#isolate significant rows for col2
	sig2 = dat[!is.na(dat[,col2]), ]
	sig2 = sig2[sig2[,col2] <= CUT,]
	
	#isolate non-significant rows for col2 (include NAs)
	notsig2 = dat[dat[,col2] > CUT,]
	
	
	#check we have all genes
	sum(nrow(sig2), nrow(notsig2))
	nrow(dat)
	
	
	#isolate rows significant in both
	sigboth = nrow(sig1[sig1[,shared.col] %in% sig2[,shared.col] , ])
	
	#isolate genes significant in only col1
	sig1only = nrow(sig1[!sig1[,shared.col] %in% sig2[,shared.col], ])
	sum(c(sigboth, sig1only))
	nrow(sig1)
	
	#isolate genes significant in only col2
	sig2only = nrow(sig2[!sig2[,shared.col] %in% sig1[,shared.col], ])
	sum(sigboth, sig2only)
	nrow(sig2)
	
	#isolate genes significant in niether
	notsig1 = dat[!dat[,shared.col] %in% sig1[,shared.col],]
	notsigBoth = nrow(notsig1[!notsig1[,shared.col] %in% sig2[,shared.col],])
	
	sum(c(notsigBoth, sig2only))
	nrow(notsig1)
	
	sum(notsigBoth, sig1only)
	nrow(notsig2)
	
	sum(sigboth, sig1only, sig2only, notsigBoth)
	
	
	sig2count = c(sigboth,sig2only)
	notsig2count = c(sig1only,notsigBoth)
	tab = as.table(rbind(sig2count, notsig2count))
	colnames(tab) = c('sig1count', 'notSig1count')
	print(tab)
	
	
	#check all genes are accounted for
	print("ARE ALL GENES ACCOUTNED FOR?")
	print("sum(tab) == nrow(dat):")
	print(sum(tab) == nrow(dataIn))
	
	#do fisher test
	ft=fisher.test(tab, alternative = 'greater')
	print(ft)
	return(ft)
}



mod.plotPCA <- function (object, intgroup = "condition", ntop = 25000, returnData = F, pcs = 10, pc1 = 1, pc2 = 2, main = "\n") 
{
    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
        length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(colData(object)))) {
        stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(object)[, intgroup, 
        drop = FALSE])
    group <- if (length(intgroup) > 1) {
        factor(apply(intgroup.df, 1, paste, collapse = " : "))
    }
    else {
        colData(object)[[intgroup]]
    }
    d <- data.frame(pca$x[,1:pcs], group = group, 
        intgroup.df, name = colnames(object))
    attr(d, "percentVar") <- percentVar[1:2]
    g = ggplot(data = d, aes_string(x = paste('PC', pc1, sep = ''),
    y = paste('PC', pc2, sep = ''), color = "group")) + 
        geom_point(size = 3) + xlab(paste0(paste0(paste0("PC", pc1), ": "), round(percentVar[pc1] * 
        100), "% variance")) + ylab(paste0(paste0(paste0("PC", pc2), ": "), round(percentVar[pc2] * 
        100), "% variance")) + coord_fixed()
   g = g + ggtitle(main)
   g = g + theme_bw()
   print(g)
   return(d)
}


printPCA = function(dat, pc1, pc2){
	ggplot(data = dat, aes_string(x = paste('PC', pc1, sep = ''), y = paste('PC', pc2, sep = ''), color = "group")) + 
        geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] * 
        100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] * 
        100), "% variance")) + coord_fixed()
}




#modified version of the function provided in the DESeq package that can be run from a dataframe
#instead of DESeqTransform object output from
mod.plotPCA.df <- function (df, coldat, intgroup = "condition", ntop = 25000, returnData = F, pcs = 10, pc1 = 1, pc2 = 2, main = "\n", SIZE = 5, SHAPE = 19) 
{
    rv <- rowVars(df)
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
        length(rv)))]
    pca <- prcomp(t(df[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    intgroup.df <- as.data.frame(coldat[, intgroup, 
        drop = FALSE])
    group <- if (length(intgroup) > 1) {
        factor(apply(intgroup.df, 1, paste, collapse = " : "))
    }
    else {
        coldat[[intgroup]]
    }
    d <- data.frame(pca$x[,1:pcs], group = group, 
        intgroup.df, name = colnames(df))
    attr(d, "percentVar") <- percentVar[1:2]
    g = ggplot(data = d, aes_string(x = paste('PC', pc1, sep = ''),
    y = paste('PC', pc2, sep = ''), color = "group")) + 
        geom_point(size=SIZE, shape=SHAPE) + xlab(paste0(paste0(paste0("PC", pc1), ": "), round(percentVar[pc1] * 
        100), "% variance")) + ylab(paste0(paste0(paste0("PC", pc2), ": "), round(percentVar[pc2] * 
        100), "% variance")) + coord_fixed()
   g = g + ggtitle(main)
   g = g + theme_bw()
   print(g)
   if (returnData == T){
   return(d)
   }
}


draw_mean_line=function(dat, color){
	mns = tapply(dat$log2FoldChange, dat$num, mean)
	xs=as.numeric(names(mns))
	pdat = data.frame(mns, xs)
	lines(mns~xs, data = pdat[pdat$xs < 30,], col = color, lwd = 2)
}


#function to output a dataframe of the trait correlations for a given module
get_module_corrs = function(m, TRAIT, ptable){
	genes = rownames(geneModuleMembership)
	column = match(m, modNames);
	moduleGenes = moduleColors==m;
	res=data.frame(abs(geneModuleMembership[moduleGenes, column]), abs(geneTraitSignificance[moduleGenes, 1]))
	colnames(res) = c('kme', 'trait_corr')
	rownames(res) = genes[moduleColors==m]
	annots= ptable[ptable$locusName %in% rownames(res), c('genbank.prot', 'protein.name', 'locusName')]
	annots=annots[!duplicated(annots$locusName),]
	rownames(annots) = annots$locusName
	res=merge(res, annots, by=0, all.x=T)
	res=res[order(res$kme, decreasing = T),]
	res$locusName<-NULL
	colnames(res)[1] <- 'locusName'
	return(res)
}



#function to set up same color scheme as ggplot2 default
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}





#plot correlations
do.plot = function(dat, y, x, ylab, xlab, xlim=NULL, ylim = NULL, main = '\n'){
	plot(dat[,y]~dat[,x],
	xlab = xlab,
	ylab = ylab,
	ylim = ylim,
	xlim = xlim,
	main = main,
	col = 'grey50')
	lm1 = lm(dat[,y]~dat[,x])
	print(summary(lm1))
	print(cor.test(dat[,y], dat[,x], method = 'spearman'))
	abline(lm1, col = 'grey')
}



#merge_protein_names
merge_protein_names = function(dat, sort.col=F, DECREASING=F){
	load('datasets/ProteinTable.Rdata')
	dat$locusName=rownames(dat)
	x=merge(ptable, dat, by='locusName')
	rownames(x)=x$locusName
	x$contig<-NULL
	x$locusName<-NULL
	if(sort.col != FALSE){
		x=x[order(x[,sort.col], decreasing= DECREASING),]
	}
	return(x)
}






#function to calculate z-score for a vector
zscore = function(x){
	mnx = mean(x)
	sdx=sd(x)
	z=(x-mnx)/sdx
	return(z)
}



get_ld_diffs = function(trait.df, ld.col, geno1, geno2){
	trait.df= trait.df[!is.na(trait.df[,ld.col]),]
	mns=tapply(trait.df[,ld.col], trait.df[,'treat'], mean)
	home=mns[geno1]
	sub2=trait.df[trait.df$treat==geno2,]
	sub2$acclim=abs(home-sub2[,ld.col])
	return(sub2)
}








get_acclim=function(trans, home){
	trans$acclim=abs(home-trans$gbm.ld1)
	return(trans)
}

YLINE.POS = 2.5

plot_acclim_fitness = function(df, fit.proxy, nclim, plot.subs=F, XLAB='Match Score', YLAB=fit.proxy, YLIM=NULL){
	LTY=1
	colors=get.cols(df$treat)
	plot(df[, fit.proxy]~df[, nclim], pch=21, bg=colors, col='black', cex=2, mgp=MGP, xlab=XLAB, ylab='', cex.axis= CEX.AXIS, cex.lab=CEX.AXIS, ylim=YLIM)
	title(ylab=YLAB, line= YLINE.POS, cex.lab=CEX.AXIS)
	df=df[!is.na(df[,fit.proxy]),]
	#plot regressions for subsets
	if (plot.subs != F){
		print("Plotting Subsets")
		sub.lwd=1
		sub.ko = df[df$treat=="KO" & !is.na(df[,fit.proxy]),]
		sub.ok = df[df$treat=="OK" & !is.na(df[,fit.proxy]),]
		lm.ko=lm(sub.ko[, fit.proxy]~ sub.ko[, nclim])
		lm.ok=lm(sub.ok[, fit.proxy]~ sub.ok[, nclim])
		clip(x1=min(sub.ko[, nclim]), max(sub.ko[, nclim]), y1=min(sub.ko[,fit.proxy]), y2=max(sub.ko[,fit.proxy]))
		print(sub.ko)
		abline(lm.ko, lty=2, col=color.set[2], lwd= sub.lwd)
		# abline(v=c(min(sub.ko[, nclim]), max(sub.ko[, nclim])), col=color.set[2])
		clip(x1=min(sub.ok[, nclim]), max(sub.ok[, nclim]), y1=min(sub.ok[,fit.proxy]), y2=max(sub.ok[,fit.proxy]))
		abline(lm.ok, lty=2, col=color.set[3], lwd= sub.lwd)
		# abline(v=c(min(sub.ok[, nclim]), max(sub.ok[, nclim])), col=color.set[3])
	}

	#plot for combination
	clip(x1=min(df[, nclim]), max(df[, nclim]), y1=min(df[,fit.proxy]), y2=max(df[,fit.proxy]))
	lm1=lm(df[, fit.proxy]~df[, nclim])
	abline(lm1, lwd=2, lty=LTY)
	p=round(summary(lm1)$coefficients[2,4], digits=3)
	R2=round(summary(lm1)$r.squared, digits=2)
	# title(main=paste('\n\np=', p))
	# title(main=paste('R2=', R2))
	title(main=bquote("R"^2 ~ "=" ~ .(R2)))
	title(main=bquote("p =" ~ .(p)), line=.75)
	
	#output stats
	if (plot.subs != F){
		r2ko=round(summary(lm.ko)$r.squared, digits=2)
		r2ok=round(summary(lm.ok)$r.squared, digits=2)
		pko=round(summary(lm.ko)$coefficients[2,4], digits=3)
		pok=round(summary(lm.ok)$coefficients[2,4], digits=3)
		R2s=c(r2ko, r2ok, R2)
		ps = c(pko, pok, p)
		group=c('KO', 'OK', 'Combined')
		print('----------------------');print('Include Origin As a Factor:')
		lm2=lm(clim.df[,fit.proxy] ~ clim.df$z + clim.df$ori)
		print(summary(lm2))
		return(data.frame(group, R2s, ps))
		
	}
	else{print(summary(lm1))}
}




plot_sub_lm = function(df, fit.proxy, ld.col, sub.treatment, color, CEX=2, LWD=3, LTY=1, PLOT=T, LINE=T, line.col=color){
	sub = na.omit(df[df$treat== sub.treatment, c(fit.proxy, ld.col)])
	lm.sub = lm(sub[,fit.proxy]~sub[, ld.col])
	if(PLOT){
		points(sub[,fit.proxy]~sub[,ld.col], bg=color, cex=CEX, pch=21)
		}
	if (LINE ==TRUE){
		clip(x1=min(sub[, ld.col]), max(sub[, ld.col]), y1=min(sub[,fit.proxy]), y2=max(sub[,fit.proxy]))
		abline(lm.sub, col= line.col, lty= LTY, lwd=LWD)
		clip(-1e8,1e8,-1e8,1e8)
	}
	r2=round(summary(lm.sub)$r.squared, 2)
	p=round(summary(lm.sub)$coefficients[2,4], 4)
	stats=c(sub.treatment, r2, p)
	return(stats)
	
}

plot_ld_fitness2= function(traits, fit.proxy, ld.col, legend.pos='topright', plot.natives = T, YLAB=fit.proxy, XLAB=ld.col, YLIM=NULL, XLIM=NULL){
	if (plot.natives){
		plot(traits[,fit.proxy]~traits[,ld.col], bg=color.set, pch=26, mgp=MGP, xlab=XLAB, ylab='', cex.axis= CEX.AXIS, cex.lab=CEX.AXIS, ylim=YLIM, xlim=XLIM)
		title(ylab=YLAB, line= YLINE.POS, cex.lab=CEX.AXIS)
		kk=plot_sub_lm(traits, fit.proxy, ld.col, "KK", color.set[1], LINE=F, CEX=.8)
		oo=plot_sub_lm(traits, fit.proxy, ld.col, "OO", color.set[4], LINE=F, CEX=.8)
	}
	else{
		sub=rbind(traits[traits$treat=="KO",], traits[traits$treat=="OK",])
		plot(sub[,fit.proxy]~ sub[,ld.col], bg=color.set, pch=26, mgp=MGP, xlab=XLAB, ylab='', cex.axis= CEX.AXIS, cex.lab=CEX.AXIS, ylim=YLIM, xlim=XLIM)
		title(ylab=YLAB, line= YLINE.POS, cex.lab=CEX.AXIS)
	}
	ko=plot_sub_lm(traits, fit.proxy, ld.col, "KO", color=color.set[2], line.col='deepskyblue3')
	ok=plot_sub_lm(traits, fit.proxy, ld.col, "OK", color=color.set[3], line.col='orange3')
	# legend(legend.pos, c('KK', 'KO', 'OK', 'OO'), pt.bg=color.set, pch=21)
	if (plot.natives){stat.res=rbind(kk, oo, ko, ok)}
	else{stat.res=rbind(ko, ok)}
	colnames(stat.res) = c('treat', 'R2', "p.value")
	return(stat.res)
}


plot_swapped_ld_fitness = function(traits, fit.proxy, ld.col, legend.pos='topright', mgp=MGP, YLAB=fit.proxy, XLAB=ld.col, YLIM=NULL){
	k.clones = swap_lds(traits, "K", "O", ld.col, fit.proxy)
	o.clones = swap_lds(traits, "O", "K", ld.col, fit.proxy)
	all.swapped = rbind(k.clones, o.clones)
	treat = substr(all.swapped$Colony.ID, start=1, stop=2)
	cols=get.cols(treat)
	plot(all.swapped[,fit.proxy]~all.swapped[,ld.col], pch=26, mgp=MGP, ylim=YLIM, cex.axis= CEX.AXIS, cex.lab=CEX.AXIS, ylab='', xlab=XLAB)
	title(ylab=YLAB, line= YLINE.POS, cex.lab=CEX.AXIS)
	# kk=plot_sub_lm(all.swapped, fit.proxy, ld.col, "KK", color.set[1], LINE=F, CEX=.8, PLOT=F)
	# oo=plot_sub_lm(all.swapped, fit.proxy, ld.col, "OO", color.set[4], LINE=F, CEX=.8, PLOT=F)
	ko=plot_sub_lm(all.swapped, fit.proxy, ld.col, "KO", color=color.set[2], line.col='deepskyblue3')
	ok=plot_sub_lm(all.swapped, fit.proxy, ld.col, "OK", color=color.set[3], line.col='orange3')
	# legend(legend.pos, c('KK', 'KO', 'OK', 'OO'), pt.bg=color.set, pch=21)
	stat.res=rbind(ko, ok)
	colnames(stat.res) = c('treat', 'R2', "p.value")
	return(stat.res)
}


#funciton to swab discriminant function values between clone-mates
swap_lds = function(traits, ori, trans, ld.col, fit.proxy){
	clones = traits[traits$ori == ori,]
	natives = traits[traits$tra == ori,]
	transplants = traits[traits$tra == trans,]
	m = merge(natives, transplants, by = 'geno')
	native.res = data.frame(m$Colony.ID.x, m[,paste(ld.col, 'y', sep='.')], m[,paste(fit.proxy, 'x', sep='.')])
	transplant.res = data.frame(m$Colony.ID.y, m[,paste(ld.col, 'x', sep='.')], m[,paste(fit.proxy, 'y', sep='.')])
	colnames(native.res) = c('Colony.ID', ld.col, fit.proxy)
	colnames(transplant.res) = c('Colony.ID', ld.col, fit.proxy)
	res=rbind(native.res, transplant.res)
	res$treat = substr(res$Colony.ID, start=1, stop=2)
	return(res)	
}


#function to return assigned colors from color.set
#from a vector of treatments
get.cols = function(treatment){
	colors=treatment
	colors[colors=='KK']<-color.set[1]
	colors[colors=='KO']<-color.set[2]
	colors[colors=='OK']<-color.set[3]
	colors[colors=='OO']<-color.set[4]
	return(colors)
}



personalized_ghosts = function(a, plot.order = c(1,4,2,3), ALPHA=.7, LWD=1, plot.mns=F){
	ldens=tapply(a$LD1, a$treat, density)
	mns=tapply(a$LD1, a$treat, mean)
	allx <- unlist(lapply(ldens, function(e) e$x))
	ally <- unlist(lapply(ldens, function(e) e$y))
	grp=as.factor(a$treat)
	levels(grp)
	MGP=c(2.1, .75, 0.0)
	plot(allx, ally, type = "n", xlab ="Discriminant function", ylab = "Density", mgp=MGP, axes=F)
	axis(1, mgp=MGP);axis(2, mgp=MGP)
	colors=get.cols(a$treat)
	# points(x=a$LD1, y=rep(0, nrow(a)), col=colors, pch="|")
	for (i in plot.order){
		polygon(c(ldens[[i]]$x, rev(ldens[[i]]$x)), c(ldens[[i]]$y, rep(0, length(ldens[[i]]$x))), col = alpha(color.set[i], ALPHA), lwd = LWD, border=color.set[i])
	}
	if (plot.mns !=F){
		abline(v=mns)
	}
}


#Function snp_ghosts()
#Basically same as personalized_ghosts()
#Plots distributions for individual coordinates
#From DAPC results for SNP data
snp_ghosts=function(a, plot.order=c(1,2), ALPHA=.7, LWD=1, plot.mns=F){
	ldens=tapply(a$LD1, a$treat, density)
	mns=tapply(a$LD1, a$treat, mean)
	allx <- unlist(lapply(ldens, function(e) e$x))
	ally <- unlist(lapply(ldens, function(e) e$y))
	grp=as.factor(a$treat)
	levels(grp)
	MGP=c(2.1, .75, 0.0)
	plot(allx, ally, type = "n", xlab ="Discriminant function", ylab = "Density", mgp=MGP, axes=F)
	axis(1, mgp=MGP);axis(2, mgp=MGP)
	colors=c(color.set[1], color.set[4])
	for (i in plot.order){
		polygon(c(ldens[[i]]$x, rev(ldens[[i]]$x)), c(ldens[[i]]$y, rep(0, length(ldens[[i]]$x))), col = alpha(colors[i], ALPHA), lwd = LWD, border=colors[i])
	}
	if (plot.mns !=F){
		abline(v=mns)
	}
}


#Funciton order_ld()
#merges individual coordinates from 
#discriminate analysis with other phenotypic data
order_ld = function(dapc){
	dapc$Colony.ID = rownames(dapc)
	head(dapc)
	t=data.frame(traits[,c('Colony.ID')])
	colnames(t) = c('Colony.ID')
	x = merge(t, dapc, by='Colony.ID', all.x=T)
	x=x[order(x$Colony.ID),]
	if(!sum(x$Colony.ID == traits$Colony.ID)==nrow(traits)){
		print("ERROR, ROWNAMES DON'T MATCH")
	}
	else{
		ld1=x$LD1
		return(ld1)
	}
}



get_shift = function(g){
	cid=rownames(g)
	g$geno=paste(substr(cid, start=1, stop=1), substr(cid, start=3, stop=5), sep='')
	ko=g[g$treat=='KO',]
	kk=g[g$treat=='KK',]
	ok=g[g$treat=='OK',]
	oo=g[g$treat=='OO',]
	mk=merge(kk, ko, by='geno')
	mo=merge(oo,ok,by='geno')
	mk$shift=abs(mk$LD1.x - mk$LD1.y)
	mo$shift=abs(mo$LD1.x - mo$LD1.y)
	m=rbind(mk, mo)
	m$shift = zscore(m$shift)
	return(m)
}




preacc_shift_lm = function(x, fit.proxy, colors=append(gg_color_hue(3), 'grey'), angle =0){
	lm1=lm(x[,fit.proxy] ~ x$z + x$shift)
	print(summary(lm1))
	lm2=lm(x[,fit.proxy] ~ x$z + x$shift + x$ori)
	print(summary(lm2))
	av=anova(lm2)["Sum Sq"]
	part.vars = av/sum(av)
	labs=c('Pre-Acclim', 'Shift', 'Origin', 'Residual')
	res=data.frame(labs, part.vars[,1])
	rownames(res) = res$labs
	res=res[c('Shift', 'Origin', 'Pre-Acclim', 'Residual'),]
	props = paste(round(res[,2], digits=2)*100, "%", sep="")
	labels = paste(res$labs, props)
	pie(x=res[,2],labels=labels, main=fit.proxy, col=colors, init.angle= angle, xpd=T)
	return(res)
}




#----------------------------------
plot_vanilla_pca = function(fit, XLIM=c(-1,1), YLIM=c(-1,1), CEX=1.5, invertX=1){
	res.pca=prcomp(fit, retx=TRUE, center=TRUE, scale=TRUE)
	names(res.pca)
	eig <- (res.pca$sdev)^2
	variance <- eig*100/sum(eig)
	variance[1]
	# Correlation between variables and principal components
	var_cor_func <- function(var.loadings, comp.sdev){
	  var.loadings*comp.sdev
	  }
	
	# Variable correlation/coordinates
	loadings <- res.pca$rotation
	sdev <- res.pca$sdev
	var.coord <- var.cor <- t(apply(loadings, 1, var_cor_func, sdev))
	head(var.coord[, 1:4])
	a <- seq(0, 2*pi, length = 100)
	
	#plot PCA vanilla style
	#gather and plot points
	ind.coord <- res.pca$x
	head(ind.coord[, 1:4])
	center <- res.pca$center
	scale<- res.pca$scale
	getdistance <- function(ind_row, center, scale){
	  return(sum(((ind_row-center)/scale)^2))
	  }
	d2 <- apply(fit, 1, getdistance, center, scale)
	cos2 <- function(ind.coord, d2){return(ind.coord^2/d2)}
	ind.cos2 <- apply(ind.coord, 2, cos2, d2)
	head(ind.cos2[, 1:4])
	c1 = scale(ind.coord[,1], center = T, scale = T)
	c1a = c1/max(abs(c1))
	c2 = scale(ind.coord[,2], center = T, scale = T)
	c2a = c2/max(abs(c2))
	
	#set up plotting variables
	colors = get.cols(substr(rownames(ind.coord), start=1, stop=2))
	pc1v = paste(paste("PC1 (", round(variance[1], digits = 0), sep = ""), "%)", sep = "")
	pc2v = paste(paste("PC2 (", round(variance[2], digits = 0), sep = ""), "%)", sep = "")
	
	#build the plot
	plot(c1a*invertX, c2a, bg=colors, axes = F, xlab = pc1v, ylab = pc2v, xlim =XLIM, ylim =YLIM, pch = 21, cex = CEX, mgp=MGP)
	axis(1);axis(2)
	#overlay the PCA stuff
	lines( cos(a), sin(a), type = 'l', col="grey", xlab = pc1v,  ylab = pc2v)
	arrows(0, 0, var.coord[, 1]*invertX, var.coord[, 2], length = 0.1, angle = 15, code = 2, lwd = 1)
	text(x=var.coord[,1]*invertX, y=var.coord[,2], labels=rownames(var.coord), cex = 1, adj=1)
	abline(h = 0, v = 0, lty = 2)
	# points(c1a, c2a, bg=colors, cex=CEX, pch=21)
}

add_column_data = function(dfx, dfy, mergeCol, col2add){
	x=merge(dfx, dfy, by = mergeCol, all.x=T)
	if( sum(x[,mergeCol] == dfx[,mergeCol])==nrow(dfx) ){
		print("Mering Successful, add this column to dfx")
	}
	else(
	print("Error"))
	return(x[,col2add])
}



######## BUILD GHOST PLOTS FOR ALL SETS + 6 MONTH ########
# pred.sup<-predict.dapc(dp,newdata=(t(vsds.6mo[degs10,]))) 
# names(pred.sup)
# pred.sup$assign
# colnames(a.vsd.supp)
# #must create another dataframe structure in order to plot these predicted values
# test<-dp
# test$ind.coord<-pred.sup$ind.scores
# test$posterior<-pred.sup$posterior
# test$assign<-pred.sup$assign
# test$grp<-as.factor(substr(colnames(vsds.6mo), start =1, stop=2)) #make sure o
# six=data.frame(test$ind.coord)
# six$treat=paste(substr(colnames(vsds.6mo), start =1, stop=2), '6mo')
# plus.six=rbind(a, six)
# #plot
# color.set=c('blue', 'dodgerblue', 'cyan', 'red', 'orange', 'firebrick')
# quartz()
# ggplot(data= plus.six, aes(LD1, fill=treat, color=treat)) + geom_density(alpha=0.6) + theme(panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +scale_color_manual(values=color.set) + scale_fill_manual(values=color.set)  #+ geom_vline(xintercept=mns, color=color.set)








#OLD VERSION OF PLOT LD FITNESS FUNCTION (FOR CROSS PLOTS)

# #function to plot correlation between LD and a fitness proxy
# #this is run in the DAPC scripts
# plot_ld_fitness=function(traits, fit.proxy, ld.col, YLAB=fit.proxy, YLIM=NULL, legend.pos='topright', plot.natives=T){
	# LWD=3
	# LTY=1
	# plot(traits[,fit.proxy]~traits[,ld.col], data=traits, ylab=YLAB, xlab="GBM Discriminant Function", mgp=MGP, ylim=YLIM)
	# #OO samples
	# if(plot.natives==TRUE){
		# clip(-1e8,1e8,-1e8,1e8)
		# d=traits[traits$ori=='O' & traits$tra=='O' & !is.na(traits[, ld.col]) & !is.na(traits[,fit.proxy]),]
		# CEX=1
		# points(d[,fit.proxy]~d[, ld.col], bg='red', pch=21, cex=CEX)  #OO samples
		# lmoo=lm(d[,fit.proxy]~d[, ld.col])
		# clip(x1=min(na.omit(d[, ld.col])), max(na.omit(d[, ld.col])), -1e8, 1e8)
		# # abline(lmoo, col = 'red', lty= LTY, lwd=LWD)
		# summary(lmoo)
		# clip(-1e8,1e8,-1e8,1e8)
		# #KK samples
		# d=traits[traits$ori=='K' & traits$tra=='K' & !is.na(traits[, ld.col]) & !is.na(traits[,fit.proxy]),]
		# points(d[,fit.proxy]~d[, ld.col], bg='blue', pch=21, cex=CEX)
		# lmkk=lm(d[,fit.proxy]~d[, ld.col])
		# clip(x1=min(na.omit(d[, ld.col])), max(na.omit(d[, ld.col])), -1e8, 1e8)
		# # abline(lmkk, col = 'blue', lty= LTY, lwd=LWD)
		# summary(lmkk)
		# clip(-1e8,1e8,-1e8,1e8)
		# legend(legend.pos, c('KK', 'KO', 'OK', 'OO'), pt.bg=color.set, pch=21)
		# head(traits)
	# }
	# #KO samples
	# CEX=1.5
	# d=traits[traits$ori=='K' & traits$tra=='O' & !is.na(traits[, ld.col]) & !is.na(traits[,fit.proxy]),]
	# points(d[,fit.proxy]~d[, ld.col], bg=color.set[2], pch=21, cex=CEX) #KO samples
	# lmko=lm(d[,fit.proxy]~d[, ld.col])
	# clip(x1=min(na.omit(d[, ld.col])), max(na.omit(d[, ld.col])), -1e8, 1e8)
	# abline(lmko, col = 'deepskyblue3', lty= LTY, lwd=LWD)
	# summary(lmko)
	# clip(-1e8,1e8,-1e8,1e8)
	# #OK samples
	# d=traits[traits$ori=='O' & traits$tra=='K' & !is.na(traits[, ld.col]) & !is.na(traits[,fit.proxy]),]
	# points(d[,fit.proxy]~d[, ld.col], bg=color.set[3], pch=21, cex=CEX)  #OK samples
	# lmok=lm(d[,fit.proxy]~d[, ld.col])
	# clip(x1=min(na.omit(d[, ld.col])), max(na.omit(d[, ld.col])), -1e8, 1e8)
	# abline(lmok, col = 'orange3', lty= LTY, lwd=LWD)
	# summary(lmok)
	# clip(-1e8,1e8,-1e8,1e8)
	# stat.list= list(summary(lmkk), summary(lmko), summary(lmok), summary(lmoo))
	# names(stat.list) = c('kk', 'ko', 'ok', 'oo')
	# rs=c(summary(lmkk)$r.squared, summary(lmko)$r.squared, summary(lmok)$r.squared, summary(lmoo)$r.squared)
	# ps=c(summary(lmkk)$coefficients[2,4], summary(lmko)$coefficients[2,4], summary(lmok)$coefficients[2,4], summary(lmoo)$coefficients[2,4])
	# statsTab=cbind(rs, ps)
	# colnames(statsTab) = c('R', 'P')
	# rownames(statsTab) = c('KK', 'KO', 'OK', 'OO')
	# print("---------------")
	# print("Stats:")
	# print(statsTab)
	# # return(stat.list)
# }

