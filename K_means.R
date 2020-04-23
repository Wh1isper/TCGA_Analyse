library(stats)
library(car)
k_means = kmeans(etSig$log2FoldChange,10)
cluster = k_means["cluster"]
#cluster_df = data.frame(gene_id = names(cluster),k_means = cluster)
#combine = merge(etSig,cluster_df,by="gene_id",all=FALSE)
combine = data.frame(etSig,k_means["cluster"])
