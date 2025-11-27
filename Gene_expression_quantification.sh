## set path
# folder
workdir=<SET_YOUR_PATH>
ref=<YOUR_PATH_TO_REFERENCE>
# file
reference=<YOUR_PATH_TO_REFERENCE/Oar_rambouillet_v1.0_genomic.fna
gtf=<YOUR_PATH_TO_REFERENCE/Oar_rambouillet_v1.0_genomic.gtf


# step1: filtering of raw data
cat name.txt|while read i;do      # you may need to change name.txt
java -jar $workdir/Trimmomatic-0.39/trimmomatic-0.39.jar \
PE -threads 20 -phred33 \
$workdir/raw_data/${i}_1.fq.gz $workdir/raw_data/${i}_2.fq.gz \
$workdir/clean/${i}_paired_clean_1.fq.gz \
$workdir/clean/${i}_unpair_clean_1.fq.gz \
$workdir/clean/${i}_paired_clean_2.fq.gz \
$workdir/clean/${i}_unpair_clean_2.fq.gz \
ILLUMINACLIP:$workdir/Trimmomatic/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
>> fitter.log \
2>&1
done

# step2: mapping to the genome
cat name.txt|while read i;do
$workdir/hisat2/hisat2-2.2.1/hisat2 -t -p 20 -x $workdir/ref/Oar_rambouillet_v1.0_genomic \
-1 $workdir/clean/${i}_paired_clean_1.fq.gz \
-2 $workdir/clean/${i}_paired_clean_2.fq.gz \
-S $workdir/mapping/${i}.sam \
>> hisat.log \
2>&1
done

# step3: Convert SAM to BAM and sort by chromosome
cat name.txt | while read i; do
$workdir/samtools/samtools view -S $workdir/sam/${i}.sam -b | \
$workdir/samtools/samtools sort -o $workdir/bam/${i}_sort.bam
done

# step4: build index
cat name.txt|while read i;do
$workdir/samtools/samtools \
index $workdir/bam/${i}_sort.bam
done

# step5:get read counts
cat name.txt|while read i
do
$workdir/featurcounts/subread-2.0.6-Linux-x86_64/bin/featureCounts \
-T 10 -p -a $workdir/ref/Oar_rambouillet_v1.0_genomic.gtf \
-o $workdir/counts/${i}.txt $workdir/bam/${i}_sort.bam \
>> counts.id.log \
2>&1
done