# step1: Normalized with the qsmooth
library(qsmooth)
setwd("<SET_YOUR_DIR>")
df <- read.delim('expression_matrix.txt',header=TRUE,sep='\t',row.names=1,check.names=F)
df <- df[apply(df, 1, function(row) max(row) >= 1), ]
dim(df)
write.table(df,"x.txt", sep = '\t')
group <- read.delim('group.txt',check.names=F)
group1 = group$y
qs_norm <- qsmooth(object = df, 
                    group_factor = group1)
write.table(qsmoothData(qs_norm),"qsmooth.txt", sep = '\t')

# step2: Sort by tissues
setwd("<SET_YOUR_DIR>")
data <- read.delim("qsmooth.txt",header=TRUE,sep='\t',row.names=1,check.names=F)
all_columns <- colnames(data)
organ_names <- c("hypothalamus", "hypophysis", "pineal","skin","mammary","rumen","abomasum","cecum","jejunum","heart",
"liver","spleen","lymph","lung","kidney","perirenal","mesentery","tail","subcutaneous", "muscle","testis","ovary","cartilage")
sorted_columns <- character()
for (organ in organ_names) {
  matching_columns <- grep(organ, all_columns, value = TRUE)
  sorted_columns <- c(sorted_columns, matching_columns)
}
sorted_data <- data[, c(sorted_columns)]
write.table(sorted_data, "sorted_qsmooth.txt", sep = '\t', col.names = TRUE, row.names = TRUE)

# step3: Calculate the average value between repeats of the same tissue type, and use TBtools to calculate the TAU value for each gene.
