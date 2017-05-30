#bisulfite_functions.R
MGP=c(2.1,1,0)
plot.mean.meth = function(dat, locus, treat, variable, MAIN){
	x=dat[dat[,'locus']==locus,]
	o=x[x[,treat]=='O',]
	k=x[x[,treat]=='K',]
	means = c(mean(k[, variable]), mean(o[, variable]))
	stdvs = c(sd(k[, variable]), sd(o[, variable]))
	sterrs=c(std.error(k[, variable]), std.error(o[, variable]))
	res=data.frame(means, stdvs, sterrs)
	rownames(res) = c("K", "O")
	xs=c(-.5, .5)
	print(means)
	print(sterrs)
	ymin=min(means) - (max(sterrs) + max(sterrs)/10)
	ymax=max(means) + (max(sterrs) + max(sterrs)/10)
	if(ymin < 0){
		ymin=0
	}
	plot(means~xs, xlim = c(-1, 1), ylim=c( ymin, ymax), axes=F, xlab="Reef", ylab = "Mean % Methylated")
	plotCI(xs, y=means, uiw=sterrs, add=T)
	axis(1, at = xs, labels= rownames(res));axis(2)
	tres = t.test(k[, variable], o[, variable])
	p = round(tres$p.value, digits=3)
	title(main=paste("\n\np=", p))
	return(res)
}



# locus = 'locus23'
# treat='trans'
# variable='methPct'
plot.mean.pair.diff = function(dat, locus, treat, variable, MAIN, higher){
	diffs = c()
	sub=dat[dat[,'locus']==locus,]
	genotypes = unique(as.character(sub$geno))
	print(genotypes)
	cumulative.diffs = c()
	for (g in genotypes){
		x = sub[sub$geno==g,]
		o=x[x[,treat]=='O',]
		k=x[x[,treat]=='K',]
		print("--------------------")
		print(g)
		mdat = merge(k, o, by = 'cpgSite')
		diffs = mdat$methPct.x - mdat$methPct.y
		print(diffs)
		cumulative.diffs = append(cumulative.diffs, diffs)
	}
	mn.diff = mean(cumulative.diffs)
	serr = std.error(cumulative.diffs)
	plotCI(x=0, y=mn.diff, uiw=serr, add = F, axes=F, ylim = c(-1,1), ylab="K clones - O clones", xlab=locus)
	axis(2)
	abline(h=0, lty=2, col='grey')
	t=t.test(cumulative.diffs, alternative='greater')
	title(main=paste('\n\np =', round(t$p.value, digits=3)))
}




gene.meth.barplot = function(dat, locus, treat, MAIN, plus=1/2){
	cols=c('dodgerblue', 'firebrick')
	x=dat[dat[,'locus']==locus,]
	# x=x[x$methPct > 0,]
	o=x[x[,treat]=='O',]
	k=x[x[,treat]=='K',]
	o$cpgSite<-factor(o$cpgSite)
	k$cpgSite<-factor(k$cpgSite)
	print("----------------------")
	print(x)
	omns=tapply(o$methPct, o$cpgSite, mean, simplify=T)
	kmns=tapply(k$methPct, k$cpgSite, mean, simplify=T)
	oerrs=tapply(o$methPct, o$cpgSite, std.error, simplify=T)
	kerrs=tapply(k$methPct, k$cpgSite, std.error, simplify=T)
	opos=tapply(o$position, o$cpgSite, mean)
	kpos=tapply(k$position, k$cpgSite, mean)
	in.both = names(omns)[names(omns) %in% names(kmns)]
	omns=omns[in.both];kmns=kmns[in.both]
	oerrs=oerrs[in.both];kerrs=kerrs[in.both]
	opos=opos[in.both];kpos=kpos[in.both]
	res=as.table(rbind(kmns, omns))
	#run t-tests
	cpgs=names(omns)
	pvalues = c()
	for (cpg in cpgs){
		subo = o[o$cpgSite == cpg, 'methPct']
		subk = k[k$cpgSite == cpg, 'methPct']
		if (length(subk) < 2){
			pvalues = append(pvalues, 1)
			next
		}
		if (length(subo) < 2){
			pvalues = append(pvalues, 1)
			next
		}
		tres = t.test(subk, subo)
		p = tres$p.value
		pvalues = append(pvalues, p)
	}
	symbol=pvalues
	symbol[symbol<0.001]<-'***'
	symbol[as.numeric(symbol)<0.01]<-'**'
	symbol[as.numeric(symbol)<0.05]<-'*'
	symbol[as.numeric(symbol)<=0.1]<-'&'
	symbol[as.numeric(symbol)>0.1]<-''
	print("--------Stats:")
	print(pvalues)
	print(symbol)
	all.errs = c(kerrs, oerrs)
	all.errs=all.errs[!is.na(all.errs)]
	maxy=max(c(kmns, omns)) + max(all.errs)
	print(kmns)
	print(omns)
	print(kerrs)
	print(oerrs[!is.na(oerrs)])
	print('----------omit----------')
	print(oerrs[!is.na(oerrs)])
	print(class(oerrs))
	print(c(as.numeric(kerrs), as.numeric(oerrs)))
	print(maxy)
	bp=barplot(res, beside=T, col = cols, ylim=c(0, maxy ), xlab=paste(locus, 'position'), ylab='Mean % Methylation', xaxt="n")
	xs = apply(bp, 2, mean) + plus
	ys = -maxy/25
	text(x=xs, y=ys, kpos, xpd=T, srt=45, pos=2)
	#legend('topright', c('K', 'O'), fill=cols)
	plotCI(x=bp[1,], y=kmns, uiw=kerrs, add=T, pch=26, xpd=T)
	plotCI(x=bp[2,], y=omns, uiw=oerrs, add=T, pch=26, xpd=T)
	sym.xs = apply(bp, 2, mean)
	heights = apply(res, 2, max)
	maxerrs = apply(rbind(oerrs, kerrs), 2, max)
	heights = heights + maxerrs
	heights = heights + .1*maxerrs
	par(xpd=T)
	text(sym.xs, heights, labels=symbol)
	par(xpd=F)
	return(list(res, kpos, bp))
}

meth_eigen_cor = function(dat, locus, module){
	moddat = data.frame( rownames(mEigs), mEigs[,c(module)])
	rownames(moddat) = rownames(mEigs)
	colnames(moddat) = c('sample' ,module)
	x=dat[dat[,'locus']==locus,]
	x$cpgSite<-factor(x$cpgSite)
	cpgs = unique(x$cpgSite)
	for (cpg in cpgs){
		sub = x[x$cpgSite==cpg,]
		mdat = merge(moddat, sub, by = 'sample')
		plot(mdat$methPct~mdat[,module], main = paste(locus, substr(cpg, start=9, stop=11), sep="\n"), xlab="Module Loading", ylab = "% Methylation", mgp=MGP)
		lm1=lm(mdat$methPct~mdat[,module])
		abline(lm1, col='red')
	}
}


mean_meth_eigen_cor = function(dat, locus, module){
	moddat = data.frame( mEigs[,c(module)])
	rownames(moddat) = rownames(mEigs)
	colnames(moddat) = c(module)
	print(head(moddat))
	x=dat[dat[,'locus']==locus,]
	sample = unique(x$sample)
	print(head(x))
	mns = data.frame(tapply(x$methPct, x$sample, mean, simplify=T))
	print(mns)
	res=merge(mns, moddat, by=0)
	print(head(res))
	colnames(res) = c('sample', 'means', 'mod.loads')
	plot(res$means~res$mod.loads, main = paste(module, locus, sep="\n"), xlab = "Module Loadings", ylab = "Mean % Methylation", mgp=MGP)
	lm1=lm(res$means~res$mod.loads)
	print(summary(lm1))
	abline(lm1, col='red')
}



emmix = function(x, means,sigma,lambda,k){
  mod <- normalmixEM(x, mu = means, sigma = sigma, lambda = lambda, k=k, arbvar=T)
  summary(mod)
  par(xpd = T)
  plot(mod, which = 2, breaks = 30, density = TRUE, title = "\n") 
  return(mod)
}


