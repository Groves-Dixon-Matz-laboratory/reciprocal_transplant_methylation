#frequentist_analyses_dapc_phyio.R



setwd("~/gitreps/reciprocal_transplant_methylation/scripts/bio_stats_course_final_project")
source('scripts/general_functions.R')
source('scripts/frequentist_functions.R')


#---------- SECTION 1:  load DAPC and physiological measures ----------#
lnames = load("datasets/dapc_trait_data.Rdata")
lnames
head(d)

#reduce to only the transplanted coral fragments
d=rbind(d[d$treat=='KO',], d[d$treat=='OK',])



#---------- SECTION 2: plot simple linear models for match-scores and physiological measures ----------#
#Select which predictor variable to use
x='snp.match'
x='gbm.match'
x='ge.match'

#plot linear models
quartz()
plot.subs=T
par(mfrow=c(2,3))
plot_acclim_fitness(d, 'GAIN', x, plot.subs=plot.subs)
plot_acclim_fitness(d, 'CARB', x, plot.subs=plot.subs)
plot_acclim_fitness(d, 'PROTEIN', x, plot.subs=plot.subs)
plot_acclim_fitness(d, 'LIPID', x, plot.subs=plot.subs)
# plot_acclim_fitness(d, 'ZOOX', x, plot.subs=plot.subs)
plot_acclim_fitness(d, 'fitpc1', x, plot.subs=plot.subs)
par(mfrow=c(1,1))



#---------- SECTION 3: compare predictive power of origin and gbm match ----------#
#There is variation by origin (or transplant site) in the fitness proxies
#Test to make sure the gbm match score is the more important
#explanatory variable than origin for each proxy
par(mfrow=c(2,3))
for (fp in fit.proxies){
	boxplot(d[,fp]~d[,'treat'], ylab=fp)
}


#subset the dataset to remove any NA values for these variables
par(mfrow=c(1,1))
dat=na.omit(d[,c('fitpc1', 'gbm.match', 'orik')])
lm2=lm(fitpc1 ~ orik + gbm.match, data=dat)
summary(lm2)
plot_lm_coeffs(lm2)



#repeat above for all fitness proxies
quartz()
par(mfrow=c(2,3))
for (fp in fit.proxies){
	lmfp = lm(d[,fp] ~ d[,'orik'] + d[,'gbm.match'])
	print("------------------")
	print(fp)
	print(summary(lmfp))
	plot_lm_coeffs(lmfp, main=fp)
}

#results are highly similar to Bayesian estimates



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
quartz()
par(mfrow=c(2,3))
for (fp in fit.proxies){
	lmfp = lm(d[,fp] ~ d[,'pre.gbm.match'] + d[,'gbm.shift'])
	print("------------------")
	print(fp)
	print(summary(lmfp))
	plot_lm_coeffs(lmfp, main=fp)
}



#now use AIC to determine whether both components of gbm-match should be included
#m0: include gbm.match (a composite of pre-match and shift)
#m1: include only gbm shift
#m2: include only pre gbm match
#m3: include gbm shift and gbm prematch




#prepare a dataset with no missing data
fit.proxy = 'fitpc1'
dat = na.omit(d[,c(fit.proxy, 'gbm.shift', 'pre.gbm.match', 'gbm.match')])


#run the models
fm0 = lm(dat[,fit.proxy] ~ dat$gbm.match)
fm1 = lm(dat[,fit.proxy] ~ dat$gbm.shift)
fm2 = lm(dat[,fit.proxy] ~ dat$pre.gbm.match)
fm3 = lm(dat[,fit.proxy] ~ dat$gbm.shift + dat$pre.gbm.match)



#compare them using AIC
model.list = list(fm0, fm1, fm2, fm3)
model.names = c('fm0', 'fm1', 'fm2', 'fm3')
f.gbm.breakdown = f.compare(model.list, model.names)
f.gbm.breakdown
par(mar=c(5, 4, 4, 2) + 0.1)
barplot(as.matrix(f.gbm.breakdown[,'weight']), beside=F, main="AIC weight")
plot.new()
legend('center', c('shift', 'pre', 'convergence', 'pre + shift'), pch=22, pt.bg=rev(gray.colors(4)), pt.cex=2)


#in agreement with bayesian methods, the top model 
#is the one that includes both shift and pre match


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

#run models
f.m00 = lm(fitpc1 ~ gbm.match, data=dat)
f.m0s = lm(fitpc1 ~ gbm.shift, data=dat)
f.m0p = lm(fitpc1 ~ pre.gbm.match, data=dat)
f.m0 = lm(fitpc1 ~ pre.gbm.match + gbm.shift, data=dat)
f.m1 = lm(fitpc1 ~ snp.match, data=dat)
f.m2 = lm(fitpc1 ~ ge.match, data=dat)
f.m3 = lm(fitpc1 ~ orik, data=dat)
f.m4 = lm(fitpc1 ~ snp.match + pre.gbm.match + gbm.shift, data=dat)
f.m5 = lm(fitpc1 ~ ge.match + pre.gbm.match + gbm.shift, data=dat)
f.m6 = lm(fitpc1 ~ orik + pre.gbm.match + gbm.shift, data=dat)
f.m7 = lm(fitpc1 ~ snp.match + orik + pre.gbm.match + gbm.shift, data=dat)
f.m8 = lm(fitpc1 ~ snp.match + ge.match + pre.gbm.match + gbm.shift, data=dat)
f.m9 = lm(fitpc1 ~ ge.match + orik + pre.gbm.match + gbm.shift, data=dat)
f.m10 = lm(fitpc1 ~ snp.match + ge.match + orik + pre.gbm.match + gbm.shift, data=dat)
f.m11 = lm(fitpc1 ~ snp.match + orik, data=dat)
f.m12 = lm(fitpc1 ~ snp.match + ge.match + orik, data=dat)

f.m6si = lm(fitpc1 ~ orik + pre.gbm.match + gbm.shift + orik*gbm.shift, data=dat)
f.m6pi = lm(fitpc1 ~ orik + pre.gbm.match + gbm.shift + orik*pre.gbm.match, data=dat)
f.m6bi = lm(fitpc1 ~ orik + pre.gbm.match + gbm.shift + orik*gbm.shift + orik*pre.gbm.match, data=dat)


summary(f.m00)
summary(f.m0s)
summary(f.m0p)
summary(f.m0)
summary(f.m1)
summary(f.m2)
summary(f.m3)
summary(f.m4)
summary(f.m5)
summary(f.m6)
summary(f.m7)
summary(f.m8)
summary(f.m9)
summary(f.m10)
summary(f.m11)
summary(f.m12)

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

#Compare the models using AIC
model.list=list(f.m1, f.m2, f.m3, f.m4, f.m5, f.m6, f.m7, f.m8, f.m9, f.m10, f.m11, f.m12)
model.names=c('f.m1', 'f.m2', 'f.m3', 'f.m4', 'f.m5', 'f.m6', 'f.m7', 'f.m8', 'f.m9', 'f.m10', 'f.m11', 'f.m12')
f.fitness.models = f.compare(model.list, model.names)
f.fitness.models



#look at how the posterior distributions for parameter values vary accross models
#to assess whether the those for gbm match are stable accross models
model.list = list(f.m0, f.m4, f.m5, f.m6, f.m7, f.m8, f.m9, f.m10)
model.names = c('f.m0', 'f.m4', 'f.m5', 'f.m6', 'f.m7', 'f.m8', 'f.m9', 'f.m10')
f.coeff = f.coeftab(model.list, model.names)
print(f.coeff$coef)
par(mfrow=c(1,1))
plot_estimates(f.coeff)


#as with the Bayesian methods, the coefficients for pre.gbm.match and gbm.shift stably > 0.



#based on these results, added m6.si, m6.pi, and m6bi with interaction terms
model.list = list(f.m00, f.m0, f.m1, f.m2, f.m3, f.m4, f.m5, f.m6, f.m7, f.m8, f.m9, f.m10, f.m6si, f.m6pi, f.m6bi)
model.names = c('f.m00', 'f.m0', 'f.m1', 'f.m2', 'f.m3', 'f.m4', 'f.m5', 'f.m6', 'f.m7', 'f.m8', 'f.m9', 'f.m10', 'f.m6si', 'f.m6pi', 'f.m6bi')
f.fitness.models = f.compare(model.list, model.names)
f.fitness.models


#based on these results, added m6.si, m6.pi, and m6bi with interaction terms
model.list = list(f.m0, f.m6pi, f.m6bi)
model.names = c('f.m0', 'f.m6pi', 'f.m6bi')
f.fitness.models = f.compare(model.list, model.names)
f.fitness.models


#based on AIC, the optimal model is m6pi 
#(model Fitness PC1 as a function of pre.gbm.match, gbm.shift, and the interaction between origin and pre.gbm.match)
#This suggests that the additional predictor variables, (gene expression match, and SNP match) provided little
#predictive power in addition to that descibed by gbm measures. The interaction factor suggests that the importance
#of already matching the target site varies by which site the fragment was transplanted to. Examining the plots
#suggests that this is the case. 

summary(f.m6pi)

par(mfrow=c(1,3))
plot_acclim_fitness(d, 'fitpc1', 'pre.gbm.match', plot.subs=T, XLAB='pre match score')
plot_acclim_fitness(d, 'fitpc1', 'gbm.shift', plot.subs=T, XLAB='pre match score')
plot_acclim_fitness(d, 'fitpc1', 'gbm.match', plot.subs=T, XLAB='overall match score')











