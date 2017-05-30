#heatmapping.R
library(pheatmap)

#Plot heatmap for tagseq
load("~/Desktop/redo/GE_3mo.Rdata")
head(vsd)
colnames(vsd) = colnames(counts)
pheatmap(cor(vsd, method='spearman'), cluster_cols = T, cluster_rows = T, clustering_distance_rows = "maximum", clustering_distance_cols = "maximum")
dim(vsd)



#Plot heatmap for GBM correlation for 3 month samples
lnames = load("~/Desktop/redo/mbd_3_and_6mo.Rdata")
meth.rld.df = vsd[,-grep("3m", colnames(vsd))]
dim(meth.rld.df)
labs = sub("_2m", "", colnames(meth.rld.df))
# labs[labs=="KK4"]<-"KK4'"
# labs[labs=="KO4"]<-"KO4'"
pheatmap(cor(meth.rld.df, method='spearman'), cluster_cols = T, cluster_rows = T, clustering_distance_rows = "maximum", clustering_distance_cols = "maximum", labels_row=labs, labels_col=labs)
dim(meth.rld.df)


#Plot heatmap for GBM correlation for 3 and 6 month samples
meth.rld.df=vsd
labs = sub("_2m", "", colnames(meth.rld.df))
labs = sub("_3m", "", labs)
# labs[labs=="KK4"]<-"KK4'"
# labs[labs=="KO4"]<-"KO4'"
pheatmap(cor(meth.rld.df, method='spearman'), cluster_cols = T, cluster_rows = T, clustering_distance_rows = "maximum", clustering_distance_cols = "maximum", labels_row=labs, labels_col=labs)
dim(meth.rld.df)



