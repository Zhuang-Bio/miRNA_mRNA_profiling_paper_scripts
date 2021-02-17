#!/bin/bash -l
#SBATCH -A snic2019-8-262
#SBATCH -p core 
#SBATCH -n 10
#SBATCH --mem-per-cpu=7000
#SBATCH -t 5-00:00:00
#SBATCH --nodes=1
#SBATCH -J wb_CW

##Download the miRBase miRNA sequences(ftp://mirbase.org/pub/mirbase/CURRENT/)
grep sapiens mature.fa |wc 
grep sapiens hairpin.fa |wc
## extract the Homo sapiens
perl -alne '{if(/^>/){if(/Homo/){$tmp=1}else{$tmp=0}};next if $tmp!=1;s/U/T/g if !/>/;print }' hairpin.fa >hairpin.human.fa
perl -alne '{if(/^>/){if(/Homo/){$tmp=1}else{$tmp=0}};next if $tmp!=1;s/U/T/g if !/>/;print }' mature.fa >mature.human.fa
perl -alne '{if(/^>/){if(/Mus/){$tmp=1}else{$tmp=0}};next if $tmp!=1;s/U/T/g if !/>/;print }' mature.fa >mature.mus.fa

##fastqc 
ls *_clean.fq.gz | while read id ; do ~/biosoft/fastqc/FastQC/fastqc $id;done

##mapping 
cd /crex/proj/snic2019-8-262/miRNA/miRNA_wb_sep/wb/clean_data
ls *_clean.fq.gz | while read id
do
gunzip ${id}_clean.fa.gz
remove_white_space_in_id.pl ${id}_clean.fa > ${id}_clean_noWhitespace.fa
mapper.pl ${id}_clean_noWhitespace.fa \
-c -i -j -k AGATCGGAAGAGCACACGTCT \
-l 18 -m -p /crex/proj/snic2019-8-262/miRNA/bowtie_index/genomeGRCh38Bowtie \
-s ${id}_reads_collapsed.fa \
-t ${id}_reads_vs_refdb.arf -v -o 8
rm ${id}_clean.fa
done

##miRDeep2
ls *_reads_vs_refdb.arf | while read id
do
miRDeep2.pl CW5_reads_collapsed.fa \
/crex/proj/snic2019-8-262/lncRNA/humanGenome/GRCh38_p13_genome_noWhiteSpace.fa \
${id}_reads_vs_refdb.arf \
/crex/proj/snic2019-8-262/miRNA/TargetScan7/mature_human_noWhiteSpace.fa \
/crex/proj/snic2019-8-262/miRNA/TargetScan7/mature_mus_noWhiteSpace.fa \
/crex/proj/snic2019-8-262/miRNA/TargetScan7/hairpin_human_noWhiteSpace.fa \
-t hsa -r ${id}_miRDeep2 2> ${id}_miRDeep2.log
done