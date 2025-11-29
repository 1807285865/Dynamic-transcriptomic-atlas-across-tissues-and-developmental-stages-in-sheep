## set path
# folder
workdir=<SET_YOUR_PATH>
ref=<YOUR_PATH_TO_REFERENCE>
# file
reference=<YOUR_PATH_TO_REFERENCE/Oar_rambouillet_v1.0_genomic.fna
gtf=<YOUR_PATH_TO_REFERENCE/Oar_rambouillet_v1.0_genomic.gtf

# step1: index
$workdir/STAR/STAR-2.7.11b/bin/Linux_x86_64/STAR \
--runThreadN 6 --runMode genomeGenerate \
--genomeDir $workdir/ref/ \
--genomeFastaFiles $workdir/ref/Oar_rambouillet_v1.0_genomic.fna \
--sjdbGTFfile $workdir/ref/Oar_rambouillet_v1.0_genomic.gtf \
--sjdbOverhang 149

# step2: align
cat name.txt|while read i;do     # you may need to change name.txt
$workdir/STAR/STAR-2.7.11b/bin/Linux_x86_64/STAR \
--runThreadN 20 --readFilesCommand zcat \
--twopassMode Basic \
--genomeDir $workdir/ref/ \
--readFilesIn $workdir/clean/${i}_paired_clean_1.fq.gz $workdir/clean/${i}_paired_clean_2.fq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix $workdir/bam/${i} \
>> $workdir/bam/hisat.log \
2>&1
done

# step3: Generate configuration files for the prep step of rMATS-turbo
mkdir -p bamConfiguration_prep

ls $workdir/bam/*.bam | sed 's/_sort.bam//g' > name.txt

cat name.txt | while read i;do
ls $workdir/bam/${i}_sort.bam > ./bamConfiguration_prep/${i}.txt
done

# step4: Generate the configuration file for the post step of rMATS-turbo. 
mkdir -p bamConfiguration_post

ls $workdir/bam/*.bam | tr '\n' ',' | sed 's/,$/\n/' > ./bamConfiguration_post/b1_1019.txt

# step5: Run rMATS-turbo with the prep and post steps separated
cat name.txt | while read i;do
python $workdir/rmats.py --gtf $workdir/ref/Oar_rambouillet_v1.0_genomic.gtf \
--tmp prep --od post_1019 \
--readLength 101 --b1 bamConfiguration_prep/${i}.txt -t paired \
--anchorLength 1 --nthread 1 --libType fr-unstranded --task prep \
--variable-read-length
done

python $workdir/rmats.py \
--gtf $workdir/ref/Oar_rambouillet_v1.0_genomic.gtf \
--tmp ./prep --od ./post_1019 \
--readLength 101 --b1 bamConfiguration_post/b1_1019.txt -t paired \
--anchorLength 1 --nthread 1 --libType fr-unstranded --task post \
--variable-read-length --statoff