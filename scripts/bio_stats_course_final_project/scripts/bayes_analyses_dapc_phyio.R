#bayes_analyses_dapc_phyio.R
library('rethinking')
setwd("~/Dixon_final_project")
source('scripts/general_functions.R')
source('scripts/bayesian_functions.R')


#---------- SECTION 1: load DAPC and physiological measures ----------#
lnames = load("datasets/dapc_trait_data.Rdata")
lnames
head(d)


#reduce to only the transplanted coral fragments
rownames(d) = d$Colony.ID
d=rbind(d[d$treat=='KO',], d[d$treat=='OK',])




#---------- SECTION 2: plot simple linear models for match-scores and physiological measures ----------#
fit.proxies = c('GAIN', 'CARB', 'PROTEIN', 'LIPID', 'fitpc1')

#Select which predictor variable to use
x='snp.match'
x='gbm.match'
x='ge.match'



#set interval for plotting
bound=0.9

#choose a set of priors
prior.list = alist(
		y ~ dnorm(mu, sigma),
		mu <- a + b*x,
		a~dnorm(0, 10),
		b~dnorm(0, 10),
		sigma~dunif(0, 10)
)


#construct and plot the bayesian single variable linear models
quartz()
par(mfrow=c(2,3))
plot_bayse_single_linear(d, xcol=x, ycol='GAIN', prior.list=prior.list, bound=bound)
plot_bayse_single_linear(d, xcol=x, ycol='CARB', prior.list=prior.list, bound=bound)
plot_bayse_single_linear(d, xcol=x, ycol='PROTEIN', prior.list=prior.list, bound=bound)
plot_bayse_single_linear(d, xcol=x, ycol='LIPID', prior.list=prior.list, bound=bound)
mx=plot_bayse_single_linear(d, xcol=x, ycol='fitpc1', prior.list=prior.list, bound=bound, MAIN='Fitness PC1')
plot.new();legend('topleft', c('KO', 'OK'), pt.bg=get.cols(c('KO', 'OK')), pch=21, cex=2)


#---------- SECTION 3: compare predictive power of origin and gbm match ----------#
#There is variation by origin (or transplant site) in the fitness proxies
#Test to make sure the gbm match score is the more important
#explanatory variable than origin for each proxy
par(mfrow=c(2,3))
for (fp in fit.proxies){
	boxplot(d[,fp]~d[,'treat'], ylab=fp)
}


#subset the dataset to remove any NA values for these variables
dat=na.omit(d[,c('fitpc1', 'gbm.match', 'orik')])



#set up prior list for model with origin and gbm match as predictors
prior.list = alist(
		fitpc1 ~ dnorm(mu, sigma),
		mu <- a + B.origin*orik + B.gbm*gbm.match,
		a~dnorm(0, 10),
		B.gbm ~dnorm(0, 1),
		B.origin ~dnorm(0,1),
		sigma~dunif(0, 10)
)


#run model
m2=map(prior.list, data=dat)


#look at results to assess
precis(m2)
par(mfrow=c(1,1))
plot(precis(m2))

#repeat above for all fitness proxies
prior.list = alist(
		y ~ dnorm(mu, sigma),
		mu <- a + b1*x1 + b2*x2,
		a~dnorm(0, 10),
		b1 ~dnorm(0, 1),
		b2 ~dnorm(0,1),
		sigma~dunif(0, 10)
)


#compare origin and gbm match score for each fitness proxy
quartz()
par(mfrow=c(2,3))
for (fp in fit.proxies){
	compare_two_predictors(d, y=fp, x1='gbm.match', x2='orik', prior.list)
}


#---------- SECTION 4: compare predictive power of shift vs pre-match ----------#
#The gbm match score can be subdivided into two components:
#the pre-match, defined as the match score for the transplanted fragments' clone-mate
#and the shift, defined as the distance along the disciminant axis the transplanted
#fragment moved away from its clone-mate toward the mean for the natives of the alternative site

#pre-match likely captures genetic variation, but could also potentially serve as a readout
#of epigenetic modificaiton due to local-scale environmental variation between colonies

#shift should capture the plastic response of gbm to transplantation


#Note that shift and pre-match are negatively associated, and that Orpheus corals shifted more
par(mfrow=c(1,1))
plot(gbm.shift~pre.gbm.match, data=d, pch=21, bg=get.cols(d$treat), cex=2)
legend('topright', c('KO', 'OK'), pt.bg=get.cols(c('KO', 'OK')), pch=21, cex=1)



#compare the explanatory power of pre-match and shift for each fitness proxy
par(mfrow=c(2,3))
for (fp in fit.proxies){
	compare_two_predictors(d, y=fp, x1='pre.gbm.match', x2='gbm.shift', prior.list)
}


#now use WAIC to determine whether both components of gbm-match should be included
#m0: include gbm.match (a composite of pre-match and shift)
#m1: include only gbm shift
#m2: include only pre gbm match
#m3: include gbm shift and gbm prematch



#prepare a dataset with no missing data
fit.proxy = 'fitpc1'
dat = na.omit(d[,c(fit.proxy, 'gbm.shift', 'pre.gbm.match', 'gbm.match')])



#pick prior values for coefficient distributions
dnorm.mn = 0
dnorm.sd = 10


#Set up prior lists for each model
m0.priors = alist(
	fitpc1 ~ dnorm(mu, sigma),
	mu <- a + B.match*gbm.match,
	a~dnorm(0, 10),
	B.match ~dnorm(dnorm.mn, dnorm.sd),
	sigma~dunif(0, 10)
) 

m1.priors = alist(
	fitpc1 ~ dnorm(mu, sigma),
	mu <- a + B.shift*gbm.shift,
	a~dnorm(0, 10),
	B.shift~dnorm(dnorm.mn, dnorm.sd),
	sigma~dunif(0, 10)
)

m2.priors = alist(
	fitpc1 ~ dnorm(mu, sigma),
	mu <- a + B.pre*pre.gbm.match,
	a~dnorm(0, 10),
	B.pre~dnorm(dnorm.mn, dnorm.sd),
	sigma~dunif(0, 10)
)

m3.priors = alist(
	fitpc1 ~ dnorm(mu, sigma),
	mu <- a + B.pre*pre.gbm.match + B.shift*gbm.shift,
	a~dnorm(0, 10),
	B.pre~dnorm(dnorm.mn, dnorm.sd),
	B.shift~dnorm(dnorm.mn, dnorm.sd),
	sigma~dunif(0, 10)
)



#run the models
m0=map(m0.priors, data=dat)
m1=map(m1.priors, data=dat)
m2=map(m2.priors, data=dat)
m3=map(m3.priors, data=dat)


#compare them using WAIC
gbm.breakdown = compare(m0, m1, m2, m3)
gbm.breakdown


#plot results
par(mfrow=c(1,3))
compare_two_predictors(d, y='fitpc1', x1='pre.gbm.match', x2='gbm.shift', prior.list)
w=c(71, 26, 3, 0)
mods=c('pre + shift', 'match', 'pre', 'shift')
t=as.matrix(w)
rownames(t) = mods
barplot(t, beside=F, main='Akaike Weight')
plot.new();legend('topleft', mods, pt.bg=grey.colors(n=4), xpd=T, pch=22, cex=1.2)


#conclude that both shift and pre-match are useful in explaining 
#variation in fitness proxies. Including these separate variables
#also sheds light on the processes that lead to relative match
#of transplanted fragments. Will use each of these rather than
#gbm.match alone in downstream models.

#--------------------------- SECTION 5: build models for model selection ---------------------------#
#genetic, epigenetic, and transriptional data are available,
#which of these are useful in predicting fitness proxies under altered
#environmental conditions?

#to answer this, compare models that include different combinations of 
#these predictor variables using information criteria

#get overview of correlations between variables
pairs( ~ snp.match + gbm.match + ge.match + pre.gbm.match + gbm.shift + orik, data=d)


#select predictor variables
xvars = c('snp.match', 'ge.match', 'gbm.match', 'pre.gbm.match', 'gbm.shift', 'orik')



#remove missing data
dat = na.omit(d[,c('fitpc1', xvars)])
dat



#MODELS:

#m00: gbm.match by itself
#m0s: include only gbm shift
#m0p: include only pre gbm match
#m0: include only gbm shift and pre gbm match
#m1: include only snp match
#m2: include only ge match
#m3: include only origin
#m4: snp-match + gbm shift + pre gbm match
#m5: ge + gbm shift + pre gbm match
#m6: origin + gbm shift + pre gbm match
#m7: snp-match + origin + gbm shift + pre gbm match
#m8: snp-match + ge + gbm shift + pre gbm match
#m9: ge + origin + gbm shift + pre gbm match
#m10: snp-match + ge + origin + gbm shift + pre gbm match
#m11: snp-match + origin
#m12: snp-match + ge + origin


#pick prior values for coefficient distributions
dnorm.mn = 0
dnorm.sd = 10


#Set up priors for each model


#m00: gbm.match by itself
m00.priors=alist(
	fitpc1 ~ dnorm(mu, sigma),
	mu <- a + B.match*gbm.match,
	a~dnorm(0, 10),
	B.match ~dnorm(dnorm.mn, dnorm.sd),
	sigma~dunif(0, 10)
) 


#m0s: include only gbm shift
m0s.priors = alist(
	fitpc1 ~ dnorm(mu, sigma),
	mu <- a + B.shift*gbm.shift,
	a~dnorm(0, 10),
	B.shift~dnorm(dnorm.mn, dnorm.sd),
	sigma~dunif(0, 10)
) 



#m0p: include only pre gbm match
m0p.priors = alist(
	fitpc1 ~ dnorm(mu, sigma),
	mu <- a + B.pre*pre.gbm.match,
	a~dnorm(0, 10),
	B.pre~dnorm(dnorm.mn, dnorm.sd),
	sigma~dunif(0, 10)
) 


#m0: include only gbm shift and pre gbm match
m0.priors = alist(
	fitpc1 ~ dnorm(mu, sigma),
	mu <- a + B.pre*pre.gbm.match + B.shift*gbm.shift,
	a~dnorm(0, 10),
	B.pre~dnorm(dnorm.mn, dnorm.sd),
	B.shift~dnorm(dnorm.mn, dnorm.sd),
	sigma~dunif(0, 10)
) 


#m1: include only snp match
m1.priors = alist(
	fitpc1 ~ dnorm(mu, sigma),
	mu <- a + B.snp*snp.match,
	a~dnorm(0, 10),
	B.snp ~dnorm(dnorm.mn, dnorm.sd),
	sigma~dunif(0, 10)
) 


#m2: include only ge match
m2.priors = alist(
	fitpc1 ~ dnorm(mu, sigma),
	mu <- a + B.ge*ge.match,
	a~dnorm(0, 10),
	B.ge ~dnorm(dnorm.mn, dnorm.sd),
	sigma~dunif(0, 10)
) 


#m3: include only origin
m3.priors = alist(
	fitpc1 ~ dnorm(mu, sigma),
	mu <- a + B.ori*orik,
	a~dnorm(0, 10),
	B.ori ~dnorm(dnorm.mn, dnorm.sd),
	sigma~dunif(0, 10)
) 


#m4: snp-match + gbm shift + pre gbm match
m4.priors = alist(
	fitpc1 ~ dnorm(mu, sigma),
	mu <- a + B.snp*snp.match + B.pre*pre.gbm.match + B.shift*gbm.shift,
	a~dnorm(0, 10),
	B.snp ~dnorm(dnorm.mn, dnorm.sd),
	B.pre~dnorm(dnorm.mn, dnorm.sd),
	B.shift~dnorm(dnorm.mn, dnorm.sd),
	sigma~dunif(0, 10)
)


#m5: ge + gbm shift + pre gbm match
m5.priors = alist(
	fitpc1 ~ dnorm(mu, sigma),
	mu <- a + B.ge*ge.match + B.pre*pre.gbm.match + B.shift*gbm.shift,
	a~dnorm(0, 10),
	B.ge ~dnorm(dnorm.mn, dnorm.sd),
	B.pre~dnorm(dnorm.mn, dnorm.sd),
	B.shift~dnorm(dnorm.mn, dnorm.sd),
	sigma~dunif(0, 10)
)


#m6: origin + gbm shift + pre gbm match
m6.priors = alist(
	fitpc1 ~ dnorm(mu, sigma),
	mu <- a + B.ori*orik + B.pre*pre.gbm.match + B.shift*gbm.shift,
	a~dnorm(0, 10),
	B.ori~dnorm(0,1),
	B.pre~dnorm(dnorm.mn, dnorm.sd),
	B.shift~dnorm(dnorm.mn, dnorm.sd),
	sigma~dunif(0, 10)
)


#m7: snp-match + origin + gbm shift + pre gbm match
m7.priors = alist(
	fitpc1 ~ dnorm(mu, sigma),
	mu <- a + B.snp*snp.match + B.ori*orik + B.pre*pre.gbm.match + B.shift*gbm.shift,
	a~dnorm(0, 10),
	B.snp ~dnorm(dnorm.mn, dnorm.sd),
	B.ori ~dnorm(dnorm.mn, dnorm.sd),
	B.pre~dnorm(dnorm.mn, dnorm.sd),
	B.shift~dnorm(dnorm.mn, dnorm.sd),
	sigma~dunif(0, 10)
)


#m8: snp-match + ge + gbm shift + pre gbm match
m8.priors = alist(
	fitpc1 ~ dnorm(mu, sigma),
	mu <- a + B.snp*snp.match + B.ge*ge.match + B.pre*pre.gbm.match + B.shift*gbm.shift,
	a~dnorm(0, 10),
	B.snp ~dnorm(dnorm.mn, dnorm.sd),
	B.ge ~dnorm(dnorm.mn, dnorm.sd),
	B.pre~dnorm(dnorm.mn, dnorm.sd),
	B.shift~dnorm(dnorm.mn, dnorm.sd),
	sigma~dunif(0, 10)
)


#m9: ge + origin + gbm shift + pre gbm match
m9.priors = alist(
	fitpc1 ~ dnorm(mu, sigma),
	mu <- a + B.ge*ge.match + B.ori*orik + B.shift*gbm.shift + B.pre*pre.gbm.match,
	a~dnorm(0, 10),
	B.ge ~dnorm(dnorm.mn, dnorm.sd),
	B.ori~dnorm(dnorm.mn, dnorm.sd),
	B.pre~dnorm(dnorm.mn, dnorm.sd),
	B.shift~dnorm(dnorm.mn, dnorm.sd),
	sigma~dunif(0, 10)
)


#m10: snp-match + ge + origin + gbm shift + pre gbm match
m10.priors = alist(
	fitpc1 ~ dnorm(mu, sigma),
	mu <- a + B.snp*snp.match + B.ge*ge.match + B.ori*orik + B.shift*gbm.shift + B.pre*pre.gbm.match,
	a~dnorm(0, 10),
	B.snp~dnorm(dnorm.mn, dnorm.sd),
	B.ge ~dnorm(dnorm.mn, dnorm.sd),
	B.ori~dnorm(dnorm.mn, dnorm.sd),
	B.pre~dnorm(dnorm.mn, dnorm.sd),
	B.shift~dnorm(dnorm.mn, dnorm.sd),
	sigma~dunif(0, 10)
)


#m11: snp-match + origin
m11.priors = alist(
	fitpc1 ~ dnorm(mu, sigma),
	mu <- a + B.snp*snp.match + B.ori*orik,
	a~dnorm(0, 10),
	B.snp~dnorm(dnorm.mn, dnorm.sd),
	B.ori~dnorm(dnorm.mn, dnorm.sd),
	sigma~dunif(0, 10)
)


#m12: snp-match + ge + origin
m12.priors = alist(
	fitpc1 ~ dnorm(mu, sigma),
	mu <- a + B.snp*snp.match + B.ge*ge.match + B.ori*orik,
	a~dnorm(0, 10),
	B.snp~dnorm(dnorm.mn, dnorm.sd),
	B.ge ~dnorm(dnorm.mn, dnorm.sd),
	B.ori~dnorm(dnorm.mn, dnorm.sd),
	sigma~dunif(0, 10)
)


#m6.si: origin + gbm shift + pre gbm match + shift*origin interaction
m6.si.priors = alist(
	fitpc1 ~ dnorm(mu, sigma),
	mu <- a + B.ori*orik + B.pre*pre.gbm.match + B.shift*gbm.shift + Bos*orik*gbm.shift,
	a~dnorm(0, 10),
	B.ori~dnorm(dnorm.mn, dnorm.sd),
	B.pre~dnorm(dnorm.mn, dnorm.sd),
	B.shift~dnorm(dnorm.mn, dnorm.sd),
	Bos~dnorm(dnorm.mn, dnorm.sd),
	sigma~dunif(0, 10)
)

#m6.pi: origin + gbm shift + pre gbm match + pre-match*origin interaction
m6.pi.priors = alist(
	fitpc1 ~ dnorm(mu, sigma),
	mu <- a + B.ori*orik + B.pre*pre.gbm.match + B.shift*gbm.shift + Bop*orik*pre.gbm.match,
	a~dnorm(0, 10),
	B.ori~dnorm(dnorm.mn, dnorm.sd),
	B.pre~dnorm(dnorm.mn, dnorm.sd),
	B.shift~dnorm(dnorm.mn, dnorm.sd),
	Bop~dnorm(dnorm.mn, dnorm.sd),
	sigma~dunif(0, 10)
)


#m6.bi: origin + gbm shift + pre gbm match + shift*origin interaction + pre-match*origin interaction
m6.bi.priors = alist(
	fitpc1 ~ dnorm(mu, sigma),
	mu <- a + B.ori*orik + B.pre*pre.gbm.match + B.shift*gbm.shift + Bos*orik*gbm.shift + Bop*orik*pre.gbm.match,
	a~dnorm(0, 10),
	B.ori~dnorm(0,1),
	B.pre~dnorm(dnorm.mn, dnorm.sd),
	B.shift~dnorm(dnorm.mn, dnorm.sd),
	Bos~dnorm(dnorm.mn, dnorm.sd),
	Bop~dnorm(dnorm.mn, dnorm.sd),
	sigma~dunif(0, 10)
)


#run models
m00=map(m00.priors, data=dat)
m0s=map(m0s.priors, data=dat)
m0p=map(m0p.priors, data=dat)
m0=map(m0.priors, data=dat)
m1=map(m1.priors, data=dat)
m2=map(m2.priors, data=dat)
m3=map(m3.priors, data=dat)
m4=map(m4.priors, data=dat)
m5=map(m5.priors, data=dat)
m6=map(m6.priors, data=dat)
m7=map(m7.priors, data=dat)
m8=map(m8.priors, data=dat)
m9=map(m9.priors, data=dat)
m10=map(m10.priors, data=dat)
m11=map(m11.priors, data=dat)
m12=map(m12.priors, data=dat)
m6.si=map(m6.si.priors, data=dat)
m6.pi=map(m6.pi.priors, data=dat)
m6.bi=map(m6.bi.priors, data=dat)

precis(m00)
precis(m0)
precis(m1)
precis(m2)
precis(m3)
precis(m4)
precis(m5)
precis(m6)
precis(m7)
precis(m8)
precis(m9)
precis(m10)
precis(m11)
precis(m12)
precis(m6.pi)
precis(m6.bi)

#m00: gbm match by itself
#m0s: include only gbm shift
#m0p: include only pre gbm match
#m0: include only gbm shift and pre gbm match
#m1: include only snp match
#m2: include only ge match
#m3: include only origin
#m4: snp-match + gbm shift + pre gbm match
#m5: ge + gbm shift + pre gbm match
#m6: origin + gbm shift + pre gbm match
#m7: snp-match + origin + gbm shift + pre gbm match
#m8: snp-match + ge + gbm shift + pre gbm match
#m9: ge + origin + gbm shift + pre gbm match
#m10: snp-match + ge + origin + gbm shift + pre gbm match
#m11: snp-match + origin
#m12: snp-match + ge + origin

#--------------------------- SECTION 6: compare models ---------------------------#

#Compare the models using WAIC
fitness.models = compare(m00, m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12)
fitness.models


#look at how the posterior distributions for parameter values vary accross models
#to assess whether the those for gbm match are stable accross models
x=coeftab(m0, m4, m5, m6, m7, m8, m9, m10)
par(mfrow=c(1,1))
plot(x)
print(x)

#note that B.pre (effect of pre-acclimation) is highly stable
#B.shift is second most stable


#based on these results, added m6.si, m6.pi, and m6bi with interaction terms
fitness.models = compare(m00, m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m6.si, m6.pi, m6.bi)
fitness.models

#replicate the WAIC estimates 10x
for (i in 1:10){
	print(i)
	fitness.models = compare(m00, m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m6.si, m6.pi, m6.bi)
	print(fitness.models)
}


#replicate the WAIC estimates 10x for just the top 3 models
for (i in 1:10){
	print(i)
	fitness.models = compare(m0, m6.pi, m6.bi)
	print(fitness.models)
}

#The inconsistency of the results indicates that the top three models, (m0pi, m0bi, and m0)
#are roughly equally likely to be the optimal model. AIC returns m6pi as the optimal model.
#The inconsistency here seems to arise from uncertainty in the WAIC estimates, which are
#larger for the interaction models than the simple m0. m6.pi is still most often has the lowest
#WAIC, based on this, along with the AIC results, and the plots below, m6pi seems to be the
#the optimal model.




source('scripts/frequentist_functions.R')
par(mfrow=c(1,3))
plot_acclim_fitness(d, 'fitpc1', 'pre.gbm.match', plot.subs=T, XLAB='pre match score')
plot_acclim_fitness(d, 'fitpc1', 'gbm.shift', plot.subs=T, XLAB='overall match score')
plot_acclim_fitness(d, 'fitpc1', 'gbm.match', plot.subs=T, XLAB='overall match score')
