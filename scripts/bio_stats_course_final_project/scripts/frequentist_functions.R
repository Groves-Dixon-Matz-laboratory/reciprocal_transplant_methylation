#frequentist_functions.R
library(plotrix)
library(reshape2)


#Funciton to plot linear model for transplanted fragments for two variables
#Arguments:
#df = input dataframe
#fit.proxy = column name for response variable
#nclim = column name for predictor variable
#plot.subs = boolean for whether to plot individual linear models for two transplant groups
#XLAB = label for x axis
plot_acclim_fitness = function(df, fit.proxy, nclim, plot.subs=F, XLAB='Matching score'){
	colors=df$treat
	colors[colors=='KO']<-color.set[2]
	colors[colors=='OK']<-color.set[3]
	plot(df[, fit.proxy]~df[, nclim], pch=21, bg=colors, col='black', cex=2, mgp=MGP, xlab=XLAB, ylab=fit.proxy)
	
	#plot regressions for subsets
	if (plot.subs != F){
		sub.lwd=1
		sub.ko = df[df$treat=="KO",]
		sub.ok = df[df$treat=="OK",]
		lm.ko=lm(sub.ko[, fit.proxy]~ sub.ko[, nclim])
		lm.ok=lm(sub.ok[, fit.proxy]~ sub.ok[, nclim])
		clip(x1=min(na.omit(sub.ko[, nclim])), max(na.omit(sub.ko[, nclim])), y1=min(na.omit(df[,fit.proxy])), y2=1e8)
		abline(lm.ko, lty=2, col=color.set[2], lwd= sub.lwd)
		clip(x1=min(na.omit(sub.ok[, nclim])), max(na.omit(sub.ok[, nclim])), y1=min(na.omit(df[,fit.proxy])), y2=1e8)
		abline(lm.ok, lty=2, col=color.set[3], lwd= sub.lwd)
	}

	#plot for combination
	clip(x1=min(na.omit(df[, nclim])), max(na.omit(df[, nclim])), y1=min(na.omit(df[,fit.proxy])), 1e8)
	lm1=lm(df[, fit.proxy]~df[, nclim])
	abline(lm1, lwd=2, lty=2)
	p=round(summary(lm1)$coefficients[2,4], digits=3)
	R2=round(summary(lm1)$r.squared, digits=2)
	title(main=paste('\n\np=', p))
	title(main=paste('R2=', R2))
	
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
		lm2=lm(df[,fit.proxy] ~ df[,nclim] + df$ori)
		print(summary(lm2))
		return(data.frame(group, R2s, ps))
		
	}
	else{print(summary(lm1))}
}


plot_lm_coeffs = function(lin.mod, main=''){
	coef=summary(lin.mod)$coefficients
	par(mar=c(5, 8, 4, 2) + 0.1)
	plotCI(x=coef[,1], y=nrow(coef):1, uiw=coef[,2], err='x', axes=F, ylim = c(0, nrow(coef)+1), xlab='Estimate', ylab='', sfrac=0, main=main)
	axis(1)
	axis(2, at = nrow(coef):1, labels=rownames(coef), las=2)
	box()
	abline(h=1:nrow(coef), lwd=1, col='grey', lty=3)
	abline(v=0, lty=1, col='grey', lwd=1)
	plotCI(x=coef[,1], y=nrow(coef):1, uiw=coef[,2], err='x', axes=F, ylim = c(0, nrow(coef)+1), xlab='Estimate', ylab='', add=T, sfrac=0, lwd=1.5)
	par(mar=c(5, 4, 4, 2) + 0.1)
}



f.compare = function(model.list, model.names){
	lapply(model.list, extractAIC)
	res=data.frame(t(data.frame(lapply(model.list, extractAIC))))
	rownames(res) = model.names
	colnames(res) = c('d.f', 'AIC')
	res=res[,c('AIC', 'd.f')]
	res$dAIC = abs(res$AIC-min(res$AIC))
	daic = res$dAIC
	den = sum(exp(-.5*daic))
	res$weight = round(exp(-.5*daic)/den, digits=3)
	res=res[order(res$AIC),]
	return(res)
}




get_coeff=function(m){
	c=summary(m)$coefficients
	return(c)
}




f.coeftab = function(model.list, model.names){
	coeffs = lapply(model.list, function(x) get_coeff(x))
	res=data.frame()
	for (i in 1:length(coeffs)){
		c=data.frame(coeffs[i])
		c$model=model.names[i]
		c$predictor = rownames(c)
		colnames(c) = c('Estimate', 'Std.Error', 't.value', 'p', 'model', 'predictor')
		res=rbind(res, c)
	}
	eres = res[,c('Estimate', 'predictor', 'model')]
	eres$Estimate = round(eres$Estimate, digits=2)
	ctab = dcast(eres, predictor~model, value.var='Estimate')
	rownames(ctab) = ctab$predictor
	ctab=ctab[,model.names]
	
	seres = res[,c('Std.Error', 'predictor', 'model')]
	seres$Std.Error = round(seres$Std.Error, digits=2)
	setab = dcast(seres, predictor~model, value.var='Std.Error')
	rownames(setab) = setab$predictor
	setab=setab[,model.names]
	final = list(model.names, ctab, setab)
	names(final) = c('model.names', 'coef', 'se')
	return(final)
}



plot_estimates = function(coeff.ob){
	cs = t(data.frame(coeff.ob$coef))
	ses = t(data.frame(coeff.ob$se))
	cs
	ses
	ys= (dim(ses)[1]*dim(ses)[2]) : 1
	l1 = c()
	for (i in colnames(cs)){
		l1=append(l1, rep(i, nrow(cs)))
	}
	l2=rep(coeff.ob$model.names, ncol(cs))
	labs=paste(l1,l2, sep='   ')
	
	par(mar=c(5, 10, 4, 2) + 0.1)
	plotCI(x=cs, y=ys, uiw=ses, err='x', axes=F, sfrac=0, lwd=1, ylab='', xlab='Estimate')
	axis(1)
	axis(2, at = ys, labels=labs, las=2)
	box()
	abline(h=ys, lwd=1, col='grey', lty=3)
	abline(v=0, lty=1, col='grey', lwd=1)
	plotCI(x=cs, y=ys, uiw=ses, err='x', axes=F, sfrac=0, lwd=1.5, add=T)
	par(mar=c(5, 4, 4, 2) + 0.1)
}

