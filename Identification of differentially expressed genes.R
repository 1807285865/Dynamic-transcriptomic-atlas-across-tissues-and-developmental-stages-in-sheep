# Identification of differentially expressed genes using DESeq2
library(DESeq2)
setwd("<SET_YOUR_DIR>")
data <- read.delim("D0_M2.txt",header=TRUE,sep='\t',row.names=1,check.names=F)
condition <- factor(c(rep("D0", 8), rep("M2", 8)), levels = c("D0","M2"))
colData <- data.frame(row.names = colnames(data), condition)
dds <- DESeqDataSetFromMatrix(data, colData, design = ~condition)
dds$condition <- relevel(dds$condition, ref = "D0")
dds <- DESeq(dds)
res = results(dds, contrast=c("condition", "D0", "M2"))
DEG_DESeq2 <- as.data.frame(res)
FC <- 2
padj <- 0.05
DEG_DESeq2$Significant <- "normal"
up <- intersect(which(DEG_DESeq2$log2FoldChange > log2(FC) ),
                which(DEG_DESeq2$padj < padj))
down <- intersect(which(DEG_DESeq2$log2FoldChange < (-log2(FC))),
                  which(DEG_DESeq2$padj < padj))
DEG_DESeq2$Significant[up] <- "up"
DEG_DESeq2$Significant[down] <- "down"
table(DEG_DESeq2$Significant)
write.table(DEG_DESeq2,"D0_M2_all_genes.xls",
            row.names = T,sep="\t",quote = F)
diff_p <- DEG_DESeq2[order(DEG_DESeq2$padj),]
diff_gene_deseq2_p <-subset(diff_p,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
write.table(diff_gene_deseq2_p,file= "D0_M2_DEGs.xls",sep = '\t',row.names = T)