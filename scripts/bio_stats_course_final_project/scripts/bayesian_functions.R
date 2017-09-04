#bayesian_functions.R




plot_bayse_single_linear = function(d, xcol, ycol, prior.list, bound=0.95, MAIN=ycol){
	prep.dat = function(d){
	dat = na.omit(d[,c(xcol,ycol)])
	colnames(dat) = c('x','y')
	return(dat)
	}
	dat=prep.dat(d)

	
	print('-------------------------------')
	print(paste("X = ", xcol))
	print(paste("Y = ", ycol))
	m1 <- map(prior.list, data=dat)
	print("Precis Output:")
	print(precis(m1, prob=bound))
	x.seq <- seq(from=min(dat[,'x']), to=max(dat[,'x']), length.out=1000)
	pred.dat <- list(x=x.seq)
	mu<-link(m1, data=pred.dat)
	mu.mean = apply(mu, 2, mean)
	mu.PI <- apply(mu, 2, PI, prob=bound)
	sim.mbd<-sim(m1, data=pred.dat)
	mbd.PI<-apply(sim.mbd, 2, PI, prob=bound)
	#plot results
	plot(y~x, dat, cex=1, axes=T,mgp=MGP, main=MAIN, pch=21, bg=get.cols(substr(rownames(dat), start=1,stop=2)), xlab=xcol,ylab=ycol)
	lines(x.seq, mu.mean, col='red')
	shade(mu.PI, x.seq)
	shade(mbd.PI, x.seq)
	points(y~x, dat, cex=2, pch=21, bg=get.cols(substr(rownames(dat), start=1,stop=2)))
	return(m1)
}


compare_two_predictors = function(d, ycol, x1, x2, prior.list){
	dat = na.omit(d[,c(ycol, x1, x2)])
	colnames(dat) = c('y', 'x1', 'x2')
	m = map(prior.list, data=dat)
	print(precis(m))
	plot(precis(m))
	title(main=ycol)
}