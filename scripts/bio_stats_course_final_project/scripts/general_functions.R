#general_functions.R
#This script has functions that are called for both the frequentist and bayesian approaches

#Some global plotting variables
MGP=c(2.1, .75, 0.0)
color.set=c('blue', 'deepskyblue', 'darkgoldenrod1', 'firebrick1')
fit.proxies = c('GAIN', 'CARB', 'PROTEIN', 'LIPID', 'fitpc1')


#Funciton to assign colors based on sample names
get.cols = function(treatment){
	colors=treatment
	colors[colors=='KK']<-color.set[1]
	colors[colors=='KO']<-color.set[2]
	colors[colors=='OK']<-color.set[3]
	colors[colors=='OO']<-color.set[4]
	return(colors)
}




get_shift = function(g){
	# cid=rownames(g)
	# g$geno=paste(substr(cid, start=1, stop=1), substr(cid, start=3, stop=5), sep='')
	ko=g[g$treat=='KO',]
	kk=g[g$treat=='KK',]
	ok=g[g$treat=='OK',]
	oo=g[g$treat=='OO',]
	mk=merge(kk, ko, by='geno')
	mo=merge(oo,ok,by='geno')
	mk$shift=abs(mk$LD1.x - mk$LD1.y)
	mo$shift=abs(mo$LD1.x - mo$LD1.y)
	m=rbind(mk, mo)
	return(m)
}
