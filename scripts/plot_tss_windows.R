#plot_tss_windows.R
setwd("/Users/grovesdixon/gitreps/reciprocal_transplant_methylation/")
lnames = load("datasets/tss_window_plotting_data.Rdata") #these data generated on TACC see MBD-seq_Data_Processing_Walkthrough (LOOK AT METHYLATION DISTRIBUTION AROUDN TSSs)
lnames

#put in kb
head(out)
out$xs=out$xs/1000
out1$x1=out1$x1/1000
out2$x2=out2$x2/1000

#plot mean alone
par(mfrow=c(1,2))
ADJ=-.2
plot(mns~xs, data = out, pch = 26, ylab = "MBD-score", xlab = "Gene Position (kb)", mgp=c(2.1,1,0), axes=F)
lines(mns~xs, data = out, lwd = 2)
lines(out$mns+out$sterrs ~ out$xs, lwd=1, col='grey')
lines(out$mns-out$sterrs ~ out$xs, lwd=1, col='grey')
axis(1);axis(2);box()
legend('top', c('mean (all genes)', 'standard error'), lwd=3, col=c('black', 'grey'), inset=c(0, -.3), xpd=T)
mtext('A', side = 3, line = 1, adj = ADJ, cex = 2, xpd=T)


#plot by subset
plot(mns~xs, data = out, pch = 26, ylim = c(-4, -.25), ylab = "MBD-score", xlab = "Gene Position (kb)", mgp=c(2.1,1,0), axes=F)
axis(1);axis(2, las=2); box()
# lines(mns~xs, data = out, lwd = 2)
lines(mns1~x1, data = out1, col = 'firebrick', lwd = 2)
lines(mns2~x2, data = out2, col = 'forestgreen', lwd = 2)
legend('top', c('strongly methlyated', 'weakly methylated'), lwd=3, col=c('forestgreen', 'firebrick'), inset=c(0, -.3), xpd=T)
mtext('B', side = 3, line = 1, adj = ADJ, cex = 2, xpd=T)