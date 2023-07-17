library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

gene_names <- read.csv("gse3678ensembl.csv")
gene_names <- gene_names$To
GO_results <- enrichGO(gene = gene_names, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "ALL")
as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 15))

png("out.png", res = 250, width = 1400, height = 1800)
print(fit)
dev.off()
