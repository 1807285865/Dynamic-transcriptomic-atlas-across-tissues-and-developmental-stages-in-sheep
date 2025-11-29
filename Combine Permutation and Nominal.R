# Combine Permutation and Nominal to get the final cis-eQTL results
setwd("<SET_YOUR_DIR>")
data = read.table("165.permutation.txt.gz", header=FALSE, stringsAsFactors=FALSE)
colnames(data) = c("pid", "nvar", "shape1", "shape2", "dummy", "pval_true_df",
"sid", "dist", "ma_samples", "ma_count", "maf", "ref_factor","npval", "effects", "slope_se", "ppval", "bpval")
data$bh = p.adjust(data$bpval, method="fdr")
write.table(data[which(data$bh <= 0.05), ], "165.permutations.fdr.txt", quote=F, row.names=F, col.names=T)
data = data[which(!is.na(data[, 16])),]
set0 = data[which(data$bh <= 0.05),] 
set1 = data[which(data$bh > 0.05),]
pthreshold = (sort(set1$bpval)[1] - sort(-1.0 * set0$bpval)[1]) / 2
data$nthresholds = qbeta(pthreshold, data$shape1, data$shape2, ncp = 0, lower.tail = TRUE, log.p = FALSE)
nom=read.table("165_NG.nominals.txt.gz",header=F,stringsAsFactors=F)
colnames(nom)<-c("pid","sid", "dst", "ma_samples", "ma_count", "maf","pval", "slope", "slope_se")
nom$threshold<-data$nthresholds[match(nom$pid,data$pid)]
write.table(nom[nom$pval<nom$threshold,],"165.nominals.sigs.txt",quote=FALSE,row.names=FALSE,col.names=T,sep="\t")