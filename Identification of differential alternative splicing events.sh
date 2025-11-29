## set path
# folder
workdir=<SET_YOUR_PATH>
ref=<YOUR_PATH_TO_REFERENCE>
# file
reference=<YOUR_PATH_TO_REFERENCE/Oar_rambouillet_v1.0_genomic.fna
gtf=<YOUR_PATH_TO_REFERENCE/Oar_rambouillet_v1.0_genomic.gtf

# step1: Generate configuration files (b1.txt and b2.txt) as input files for rMATS-turbo. These two files contain comma-separated lists of BAM files for sample groups 1 and 2, respectively.
ls $prefix_dir_group1 | tr '\n' ',' | sed 's/,$/\n/' > $workdir/b1.txt 
ls $prefix_dir_group2 | tr '\n' ',' | sed 's/,$/\n/' > $workdir/b2.txt

# step2: Run rMATS-turbo to identify differential alternative splicing events.
python $workdir/rmats.py \
--b1 $workdir/b1.txt \
--b2 $workdir/b2.txt \
--gtf $workdir/ref/Oar_rambouillet_v1.0_genomic.gtf \
-t paired --nthread 20 --readLength 150 \
--od $workdir/b1_b2 \
--tmp $workdir/b1_b2

# step3: Select statistically significant events from the two-group differential alternative splicing analysis.
# Analyze in R
setwd("$workdir/b1_b2")
data <- read.table("A3SS.MATS.JC.txt", header = TRUE, sep = "\t")
filtered_data <- data[!grepl("^chrNW_", data$chr) & data$FDR < 0.05 & abs(data$IncLevelDifference) > 0.1, ]
filtered_data$regulated <- ifelse(filtered_data$IncLevelDifference > 0, "up", "down")
write.table(filtered_data, "sig_A3SS.txt", quote = FALSE, row.names = FALSE, sep = "\t")

setwd("$workdir/b1_b2")
data <- read.table("A5SS.MATS.JC.txt", header = TRUE, sep = "\t")
filtered_data <- data[!grepl("^chrNW_", data$chr) & data$FDR < 0.05 & abs(data$IncLevelDifference) > 0.1, ]
filtered_data$regulated <- ifelse(filtered_data$IncLevelDifference > 0, "up", "down")
write.table(filtered_data, "sig_A5SS.txt", quote = FALSE, row.names = FALSE, sep = "\t")

setwd("$workdir/b1_b2")
data <- read.table("RI.MATS.JC.txt", header = TRUE, sep = "\t")
filtered_data <- data[!grepl("^chrNW_", data$chr) & data$FDR < 0.05 & abs(data$IncLevelDifference) > 0.1, ]
filtered_data$regulated <- ifelse(filtered_data$IncLevelDifference > 0, "up", "down")
write.table(filtered_data, "sig_RI.txt", quote = FALSE, row.names = FALSE, sep = "\t")

setwd("$workdir/b1_b2")
data <- read.table("SE.MATS.JC.txt", header = TRUE, sep = "\t")
filtered_data <- data[!grepl("^chrNW_", data$chr) & data$FDR < 0.05 & abs(data$IncLevelDifference) > 0.1, ]
filtered_data$regulated <- ifelse(filtered_data$IncLevelDifference > 0, "up", "down")
write.table(filtered_data, "sig_SE.txt", quote = FALSE, row.names = FALSE, sep = "\t")

# step4: Run rmats2sashimiplot to generate sashimi plots for visualizing alternative splicing events of interest.
rmats2sashimiplot --b1 $bam1,$bam2,$bam3 --b2 $bam4,$bam5,$bam6 \
--event-type SE -e sashimi_events.txt --l1 b1 --l2 b2 \
--exon_s 1 --intron_s 5 -o ./output --group-info sashimi_groupInfo.txt