#nanodrop_vs_picogreen.R
setwd("/Users/grovesdixon/git_Repositories/reciprocal_transplant_methylation/bisulfite_validation/bisulfite_labwork_files")
x = read.csv('final_pico_green_results.csv')
x


plot(x$ng.ul_nanodrop~x$ng.ul, xlab = 'ng/ul picogreen', ylab = 'ng/ul nanodrop')
lm1=lm(x$ng.ul_nanodrop~x$ng.ul)
summary(lm1)
abline(lm1, col='red')


#DON'T TRUST NANO DROP!




plot(density(x[,8]))
hist(x[,8], xlab='Molarity (nM)', main='Molarity')
mol = x[,8]
summary(mol)
sd(mol)