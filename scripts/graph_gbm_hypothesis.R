#graph_gbm_hypothesis.R
#goal is to generate a graphical representation of the hypothesis that
#gene body methylation slowly tracks gene expression, and stabilizes elevated
#transcriptional state. Transcription begins with one mean, undergoes a transition
#in which mean transcription increases, then slowly decreases in variance as gene
#body methylation increases and stabilizes the new elevated transcripotional state.



#set up plotting variables
meth.sd = 0.01 #noise level in methylation
n1 = 1000    #length of phase 1
n2 = 3*n1    
transx = n1:(n1 + n1/5)

transcription.col = 'grey'
meth.col = 'black'








#-------- method 1
low.mn = 1
hi.mn = 2
sd0=0.4
sd1=0.1
sdm=0.001
xs=1:200
step=length(xs)/4


#set up time values
x1=seq(1, step)
x2=seq(step+1, step*2)
x3=seq(step*2+1, step*3)
x4=seq(step*3+1, step*4)

#set up expression values
e1 = rnorm(x1, mean=low.mn, sd = sd0)
e2 = rnorm(x2, mean=seq(low.mn, hi.mn, length.out=c(length(x2))), sd=sd0)
e3 = rnorm(x3, mean=hi.mn, sd=seq(from=sd0, to=sd1, length.out=length(x3)))
e4 = rnorm(x4, mean=hi.mn, sd=sd1)

#set up methylation values
m1 = rnorm(x1, mean=low.mn, sd = sdm)
m2 = rnorm(x2, mean=seq(low.mn, hi.mn, length.out=c(length(x2))), sd=sdm)
m3 = rnorm(x3, mean=hi.mn, sd=sdm)
m4 = rnorm(x4, mean=hi.mn, sd=sdm)

#consolidate
alle=c(e1, e2, e3, e4)
allm=c(m1, m2, m3, m4)
allx=c(x1,x2, x3, x4)
loessm = predict(loess(allm~allx, span=0.1), allx)

#plot
plot(alle~allx, pch=26, axes=F, mpg=)
lines(loessm~allx, lwd=2, col= meth.col)
lines(alle~allx, lwd=2, col= transcription.col)

#-------- method 2
low.mn = 1
hi.mn = 1.75
sd0=0.4
sd1=0.08
sdm=0.001
xs=1:200
step=length(xs)/4


#set up time values
x1=seq(1, step)
x2=seq(step+1, step*2)
x3=seq(step*2+1, step*3)
x4=seq(step*3+1, step*4)

#set up expression values
e1 = rnorm(x1, mean=low.mn, sd = sd0)
e2 = rnorm(x2, mean=hi.mn, sd=sd0)
e3 = rnorm(x3, mean=hi.mn, sd=seq(from=sd0, to=sd1, length.out=length(x3)))
e4 = rnorm(x4, mean=hi.mn, sd=sd1)

#set up methylation values
m1 = rnorm(x1, mean=low.mn, sd = sdm)
m2 = rnorm(x2, mean=seq(low.mn, hi.mn, length.out=c(length(x2))), sd=sdm)
m3 = rnorm(x3, mean=hi.mn, sd=sdm)
m4 = rnorm(x4, mean=hi.mn, sd=sdm)

#consolidate
alle=c(e1, e2, e3, e4)
allm=c(m1, m2, m3, m4)
allx=c(x1,x2, x3, x4)
loessm = predict(loess(allm~allx, span=0.1), allx)

#plot
par(mfrow=c(2,1))
plot(alle~allx, pch=26, axes=F, ylab='Transcription', xlab='');axis(2, labels=F)
lines(alle~allx, lwd=2, col= transcription.col)
plot(alle~allx, pch=26, axes=F, xlab='Time', ylab='gbM');axis(1, labels=F);axis(2, labels=F)
lines(loessm~allx, lwd=2, col= meth.col)


















#SIMULATE TRANSCRIPTION
#set up initial transcription phase where mean transcription and variance are stable
sd0 = 0.25
e1 = rnorm(n1, mean=low.mn, sd = sd0)
x1 = 1:n1
plot(e1~x1, ylim = c(low.mn-low.mn,hi.mn+3), xlim=c(0, n2), pch =26)
lines(e1~x1, ylim = c(low.mn-low.mn,hi.mn+1), col = transcription.col)

#set up transitional phase coordinates
step.size = (hi.mn-low.mn)/length(transx)
transy = seq(from=low.mn+step.size, to=hi.mn, by=step.size)


#set of second transcription phase where mean is constant, but variance
#decreases as methylation level increases
x2 = seq(max(transx)+1, n2, by=1)
TO=.05
FROM=1
step.size = (FROM-TO)/length(x2)
e2=rnorm(n2-max(transx), mean = hi.mn, sd = seq(from=FROM-step.size, to=TO, by = -step.size ))



#plot the simulated transcription data
plot(e1~x1, ylim = c(min(e1), max(e2)), xlim=c(0, n2), pch=26, axes = F, xlab = 'time', ylab='relative transcription')
lines(e1~x1, ylim = c(low.mn-low.mn,hi.mn+1), col = transcription.col)
lines(jitter(transy, amount=1.5)~transx, col= transcription.col)
lines(e2~x2, col=transcription.col)


#SIMULATE METHYLATION
#first stable phase
m1 = rnorm(n1, mean=low.mn, sd = meth.sd)

#plot the response phase where methylation slowly tracks expression
mx.all = c(transx, x2)
mx2 = mx.all[1:(length(mx.all)-length(mx.all)*.1)]
mx2.len = length(mx2)
step.size = (hi.mn-low.mn)/mx2.len
m.t2 = seq(from=low.mn+step.size, to=hi.mn, by = step.size)

#set up final stable meth stage where it has stabilized
fx = mx.all[!mx.all %in% mx2]
fm = rep(hi.mn, times=length(fx))
fmj = rnorm(length(fx), mean=hi.mn, sd = meth.sd)


all.e = c(e1, transy, e2)
all.m = c(m1, m.t2, fmj)
all.x = c(x1, transx, x2)

lines(all.e~all.x)
lines(all.m~all.x)
loessm = loess(all.m~all.x, span=0.1)
lines(all.x, predict(loessm, all.x))


#BUILD PLOT
#plot the simulated transcription data
plot(e1~x1, ylim = c(min(e2), max(e2)), xlim=c(0, n2), pch=26, axes = F, xlab = 'time', ylab='relative transcription/GBM', mgp=c(2.1,1,0))
lines(e1~x1, ylim = c(low.mn-low.mn,hi.mn+1), col = transcription.col)
lines(jitter(transy, amount=2.5)~transx, col= transcription.col)
lines(e2~x2, col=transcription.col)
#plot the simulated methylation data
loessm = loess(all.m~all.x, span=0.1)
lines(all.x, predict(loessm, all.x))


lines(m1~x1)
lines(jitter(m.t2, amount=0.02)~mx2, col= meth.col)
lines(fmj~fx)

#plot axes, legend and transition indicator
axis(1, labels=F)
axis(2, labels=F)
legend('topleft', c('transcription', 'GBM'), col=c(transcription.col, meth.col), lty=c(1,1), lwd=2)
segments(x0=min(mx2), y0=min(e2)+(min(e2)*2)+.5, y1=low.mn-2, lwd=2, lty=2, col='red')
text(x=min(mx2), y=low.mn/2, labels = 'Environmental Change', adj=-.1)





