#plot_tss_windows.R
setwd("/Users/grovesdixon/git_Repositories/reciprocal_transplant_methylation/")
lnames = load("datasets/tss_window_plotting_data.Rdata") #these data generated on TACC see MBD-seq_Data_Processing_Walkthrough (LOOK AT METHYLATION DISTRIBUTION AROUDN TSSs)
lnames

#put in kb
head(out)
out$xs=out$xs/1000
out1$x1=out1$x1/1000
out2$x2=out2$x2/1000


#plot
plot(mns~xs, data = out, pch = 26, ylim = c(-4, -.25), ylab = "Methylation", xlab = "Gene Position (kb)", mgp=c(2.1,1,0), axes=F)
axis(1);axis(2, las=2); box()
lines(mns~xs, data = out, lwd = 2)
lines(mns1~x1, data = out1, col = 'firebrick', lwd = 2)
lines(mns2~x2, data = out2, col = 'forestgreen', lwd = 2)