library(DESeq2)
library(DEFormats)


#counts <- read.delim("/Users/rishav_raj/Desktop/bioinformatics/geo/GSE_raw.tsv", sep = "\t", header = TRUE)
#rownames(counts)<-counts[,1]
#counts<-counts[,-1]
#counts

library(GEOquery)
geo_data <- getGEO("GSE58545")
annot <- geo_data[["GSE3678_series_matrix.txt.gz"]]@featureData@data
counts <- exprs(geo_data[[1]])
counts <- na.omit(counts)
counts <- round(counts)
dim(counts)
counts1 <- counts[1:25000,]
counts2 <- counts[25001:54675,]
group = rep(c("A", "B"), each = 22)
group = c(group,"B")

dge = DGEList(counts,group=group)
dds = as.DESeqDataSet(dge)

dds
dds <- DESeq(dds)
res <- results(dds)
res
head(results(dds, tidy=TRUE))
summary(res)

resOrdered <- res[order(res$padj),]
head(resOrdered)
dim(resOrdered)

sum(res$padj < 0.01, na.rm=TRUE)
plotMA(res, main="DESeq2", ylim=c(-2,2))
res.1 <- subset(resOrdered, padj < 0.1)
diffexpgenes <- sort(as.factor(rownames(res.1)))

#expressed genes with log2FoldChange>1 and pvalue<0.01
or<- res[order(res$log2FoldChange,res$pvalue),]
dim(or)
res.2 <- subset(or,log2FoldChange>1,pvalue<0.05)
dim(res.2)
expressed_genes<- sort(as.factor(rownames(res.2)))
as.data.frame(expressed_genes)


#extracting the genes with log2FC>1 and pvalue<0.01
deg<- res[which(res$log2FoldChange>1),]
deg<- deg[which(deg$pvalue<0.01 & deg$pvalue>0),]
deg<- as.data.frame(deg[,])
deg<- deg[5]
#saving the dataframe in the local
setwd("/Users/rishav_raj/Desktop/bioinformatics/geo/DEG")

# Ignore row names/numbers
write.table(deg,"de_genes3678.csv",sep=",", row.names=TRUE)


#par(mfrow=c(1,1))
#with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
#with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
#with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#plotting according to the paper
par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
with(subset(res, pvalue<0.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#PCA
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata) #using the DESEQ2 plotPCA fxn we can
