results = read.table("tables/POLYAvsTOTAL.complete.txt", header=TRUE)
rownames(results) = results$Id
head(results)


layout(matrix(1:2, 1), respect=TRUE)
plot(results$log2FoldChange, -log10(results$padj),   main="Volcano plot POLYA vs. TOTAL (ref.)", pch=".")
plot(results$log2FoldChange, -log10(results$pvalue), main="Volcano plot POLYA vs. TOTAL (ref.)", pch=".")


layout(1, respect=TRUE)
plot(results$log2FoldChange, -log10(results$pvalue),   main="Volcano plot POLYA vs. TOTAL (ref.)", pch=".", col="grey")
idx = rownames(results)[grep("HIST",rownames(results))]
text(results[idx,]$log2FoldChange, -log10(results[idx,]$pvalue), idx, cex=0.5)
points(results[idx,]$log2FoldChange, -log10(results[idx,]$pvalue), col=2)
