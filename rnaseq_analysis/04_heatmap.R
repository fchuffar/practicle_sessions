idx_genes = rownames(results)[grep("HIST",rownames(results))]
idx_samples = colnames(results)[grep("norm",colnames(results))]

# GO!
# log2 of expression for data
data = log2(as.matrix(results[idx_genes, idx_samples])+1)
# normalization... centered and reduced
data = data - apply(data, 1, mean)
data = data / apply(data, 1, sd)    
# remove rows with no variation (needed to clustering rows according to cor)
data = data[apply(data, 1, function(l) {length(unique(l))})>1, ]
# clustering samples...
# ... based on eucl. dist.
tmp_d = t(data)
d = dist(tmp_d)
hc_col = hclust(d, method="complete")
Colv = as.dendrogram(hc_col)
# clustering genes...
# ... bases on correlation
tmp_d = data
tmp_d = tmp_d[!apply(is.na(tmp_d), 1, any), ]
d = dist(1 - cor(t(tmp_d), method="pe"))
hc_row = hclust(d, method="complete")
Rowv = as.dendrogram(hc_row)      

# colors
colors=
cols = 
gplots::heatmap.2(data, 
  Rowv=Rowv, Colv=Colv, dendrogram=dendrogram, 
  col=colorRampPalette(c("cyan", "black", "red"))(20), 
  main=paste0("Mean of mean (", nrow(data), " features x ", ncol(data), " tissues)"), 
  trace="none", mar=c(10,5), useRaster=TRUE, cexRow=0.4 , cexCol=0.4
)
