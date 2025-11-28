## set path
# folder
workdir=<SET_YOUR_PATH>
ref=<YOUR_PATH_TO_REFERENCE>
# file
reference=<YOUR_PATH_TO_REFERENCE/Oar_rambouillet_v1.0_genomic.fna
gtf=<YOUR_PATH_TO_REFERENCE/Oar_rambouillet_v1.0_genomic.gtf


# step1: fastp
cat name.txt|while read i;do      # you may need to change name.txt
$workdir/fastp/fastp \
-i $workdir/raw_data/${i}_R1.fastq.gz -I $workdir/raw_data/${i}_R2.fastq.gz \
-o $workdir/cleandata/${i}_1.clean.fq.gz \
-O $workdir/cleandata/${i}_2.clean.fq.gz \
-h $workdir/cleandata/${i}_fastp_report.html \
-q 15 -l 18 -w 16
done

# step2: index
$workdir/bwa/bwa index -a bwtsw Oar_rambouillet_v1.0_genomic.fna Oar_rambouillet_v1.0

# step3: align
cat name.txt|while read i;do
$workdir/bwa/bwa \
mem -t 30 $workdir/ref/Oar_rambouillet_v1.0_genomic.fna \
$workdir/cleandata/${i}_1.clean.fq.gz \
$workdir/cleandata/${i}_2.clean.fq.gz \
| samtools sort -O bam -@ 12 -o - > $workdir/bam_raw/${i}.raw.bam
done

# step4: filtration & sort
cat name.txt|while read i;do
$workdir/samtools/samtools-1.18/samtools \
view -h -F 1804 -f 2 -q 30 $workdir/bam_raw/${i}.raw.bam | \
$workdir/samtools/samtools-1.18/samtools \
sort -n $workdir/dev/stdin -o $workdir/bam_fitter/${i}.sorted.bam
done

cat name.txt|while read i;do
$workdir/samtools/samtools-1.18/samtools \
fixmate -r $workdir/bam_fitter/${i}.sorted.bam \
$workdir/bam_fitter/${i}.fixmate.bam
done

# step5: remove PCR duplicate
cat name.txt|while read i;do
$workdir/sambamba/sambamba-1.0.1-linux-amd64-static \
markdup --remove-duplicates --nthreads 4 \
--hash-table-size=17592186044416 --overflow-list-size=20000000 --io-buffer-size=256 \
$workdir/bam_fitter/${i}.fixmate.bam \
$workdir/bam_fitter/${i}.dedup.bam
done

# step6: remove mitochondrion
cat name.txt|while read i;do
$workdir/samtools/samtools-1.18/samtools \
view --threads 4 -h $workdir/bam_fitter/${i}.dedup.bam | awk '$3!="chrM"' | \
$workdir/samtools/samtools-1.18/samtools \
sort --threads 4 /dev/stdin -o $workdir/bam_fitter/${i}.final.bam
done

# step7: index
cat name.txt|while read i;do
$workdir/samtools/samtools-1.18/samtools \
index $workdir/bam_fitter/${i}.final.bam
done

# step8: bam to bed
cat name.txt|while read i;do
bedtools bamtobed -i $workdir/bam_fitter/${i}.final.bam >  $workdir/finalbed/${i}.last.bed
done

# step9: peak calling
cat name.txt|while read i;do
macs2 callpeak -f BED -t $workdir/finalbed/${i}.last.bed \
-n $workdir/peak/${i}. \
-g 2869531331 --shift -75 --extsize 150 --nomodel -B --SPMR --q 0.01 --keep-dup all --call-summits
done
