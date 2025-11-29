## set path
# folder
workdir=<SET_YOUR_PATH>

# step1: Permutation
for j in $(seq 1 50); 
do
$workdir/bin/fastQTL.static \
--vcf $workdir/finally_output.vcf.gz \
--bed $workdir/expression_ok.bed.gz \
--cov $workdir/cov.txt.gz \
--permute 1000 10000 --normal \
--out $workdir/result_Permutation/qtl.chunk${j}.txt.gz \
--chunk $j 50&
done

zcat qtl.chunk*.txt.gz | gzip -c > 165.permutation.txt.gz

# step2: Nominal
for j in $(seq 1 50); 
do
$workdir/bin/fastQTL.static \
--vcf $workdir/finally_output.vcf.gz \
--bed $workdir/expression_ok.bed.gz \
--cov $workdir/cov.txt.gz \
--normal --chunk $j 50 \
--out $workdir/result_normal/qtl.chunk${j}.txt.gz
done

zcat qtl.chunk*.txt.gz | gzip -c > 165_NG.nominals.txt.gz

# step3: Combine Permutation and Nominal to get the final cis-eQTL results
script: "Combine Permutation and Nominal.R"

