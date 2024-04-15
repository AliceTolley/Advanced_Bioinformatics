############### STEP 1: FASTQC ###################################################################################################
#!/bin/bash

#Move to directory containing untrimmed sequences
cd ~/ngs_course/m2305732_assess/data/untrimmed_fastq

#Perform fastqc on untrimmed sequences
fastqc -t 4 *.fastq.gz

#Move fastqc files into a new directory under results
mkdir ~/ngs_course/m2305732_assess/results/fastqc_untrimmed_reads

mv *fastqc* ~/ngs_course/m2305732_assess/results/fastqc_untrimmed_reads/

#Move into results directory to unzip fastqc files
cd ~/ngs_course/m2305732_assess/results/fastqc_untrimmed_reads/

for zip in *.zip
do
unzip $zip
done

#Save fastqc summary in a text file within logs directory
cat */summary.txt > ~/ngs_course/m2305732_assess/logs/fastqc_untrimmed_summaries.txt

#Remove now redundant zipped files
rm NGS0001.R1_fastqc.zip
rm NGS0001.R2_fastqc.zip

##################################################################################################################################

####### STEP 2: TRIMMING #########################################################################################################
#Move to directory containing untrimmed sequences
cd ~/ngs_course/m2305732_assess/data/untrimmed_fastq

#Run trimmomatic to remove low quality bases with a phred score less than 33.
trimmomatic PE -threads 4 -phred33 \
 /home/ubuntu/ngs_course/m2305732_assess/data/untrimmed_fastq/NGS0001.R1.fastq.gz \
 /home/ubuntu/ngs_course/m2305732_assess/data/untrimmed_fastq/NGS0001.R2.fastq.gz \
 -baseout /home/ubuntu/ngs_course/m2305732_assess/data/trimmed_fastq/NGS0001_trimmed_R \
 ILLUMINACLIP:/home/ubuntu/anaconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:30:10 \
 TRAILING:25 MINLEN:50

#################################################################################################################################

######## STEP 3: FASTQC TRIMMED DATA ############################################################################################
#Move to directory containing trimmed sequences
cd ~/ngs_course/m2305732_assess/data/trimmed_fastq

#Perform fastqc on trimmed sequences
fastqc -t 4 NGS0001*

#Move fastqc files into new directory under results
mkdir ~/ngs_course/m2305732_assess/results/fastqc_trimmed_reads

mv *fastqc* ~/ngs_course/m2305732_assess/results/fastqc_trimmed_reads/

#Move into results directory to unzip fastqc files
cd ~/ngs_course/m2305732_assess/results/fastqc_trimmed_reads/

for zip in *.zip
do
unzip $zip
done

#Save fastqc summary in a text file within logs directory
cat */summary.txt > ~/ngs_course/m2305732_assess/logs/fastqc_trimmed_summaries.txt

#Remove now redundant zipped files
rm *.zip

#Delete untrimmed sequences to save space
rm ~/ngs_course/m2305732_assess/data/untrimmed_fastq/NGS0001.R1.fastq.gz
rm ~/ngs_course/m2305732_assess/data/untrimmed_fastq/NGS0001.R2.fastq.gz

#################################################################################################################################

########## STEP 4: ALIGNMENT DUPLICATE MARKING ##################################################################################
##Set up alignment
#Move into directory containing reference file
cd ~/ngs_course/m2305732_assess/data/reference

#Build reference index
bwa index ~/ngs_course/m2305732_assess/data/reference/hg19.fa.gz

#Make aligned_data directory
mkdir ~/ngs_course/m2305732_assess/data/aligned_data

#Run alignment on 4 threads "-t 4" to speed up alignment and complete read group information "-R"
bwa mem -t 4 -v 1 -R \
 '@RG\tID:11V6WR1.111.D1375ACXX.1.NGS0001\tSM:NGS0001\tPL:ILLUMINA\tLB:nextera-NGS0001-blood\tDT:2017-02-23\tPU:11V6WR1' \
 -I 250,50  ~/ngs_course/m2305732_assess/data/reference/hg19.fa.gz \
 ~/ngs_course/m2305732_assess/data/trimmed_fastq/NGS0001_trimmed_R_1P \
 ~/ngs_course/m2305732_assess/data/trimmed_fastq/NGS0001_trimmed_R_2P > ~/ngs_course/m2305732_assess/data/aligned_data/NGS0001.sam

#Delete trimmed sequences to save space
rm ~/ngs_course/m2305732_assess/data/trimmed_fastq/NGS0001_trimmed_R_1P
rm ~/ngs_course/m2305732_assess/data/trimmed_fastq/NGS0001_trimmed_R_1U
rm ~/ngs_course/m2305732_assess/data/trimmed_fastq/NGS0001_trimmed_R_2P
rm ~/ngs_course/m2305732_assess/data/trimmed_fastq/NGS0001_trimmed_R_2U

##Convert, sort and index SAM file
#Move into aligned_data sub-directory
cd ~/ngs_course/m2305732_assess/data/aligned_data

#Convert SAM to BAM
samtools view -h -b NGS0001.sam > NGS0001.bam

#Sort BAM file
samtools sort NGS0001.bam > NGS0001_sorted.bam

#Generate index from sorted BAM file
samtools index NGS0001_sorted.bam

##Markduplicates
picard MarkDuplicates I=NGS0001_sorted.bam O=NGS0001_sorted_marked.bam M=marked_dup_metrics.txt

samtools index NGS0001_sorted_marked.bam

##Filter SAM OR BAM
samtools view -F 1796  -q 20 -o NGS0001_sorted_filtered.bam NGS0001_sorted_marked.bam

samtools index NGS0001_sorted_filtered.bam

##Flagstat
samtools flagstat NGS0001_sorted_filtered.bam

##IdxStats
samtools idxstats NGS0001_sorted_filtered.bam

##Depth of Coverage
bedtools coverage -a NGS0001_sorted_filtered.bam -b ~/ngs_course/m2305732_assess/data/annotation.bed

##CollectInsertSizeMetrics
picard CollectInsertSizeMetrics I=NGS0001_sorted_filtered.bam O=insert_size_metrics.txt H=insert_size_metrics.pdf M=0.5
#################################################################################################################################

############ STEP 5: Variant Calling & Filtering ################################################################################
#Unzip reference file
zcat ~/ngs_course/m2305732_assess/data/reference/hg19.fa.gz > ~/ngs_course/m2305732_assess/data/reference/hg19.fa

#Build samtools index of reference
samtools faidx ~/ngs_course/m2305732_assess/data/reference/hg19.fa

#Call variants by running Freebayes
freebayes --bam ~/ngs_course/m2305732_assess/data/aligned_data/NGS0001_sorted_filtered.bam \
 --fasta-reference ~/ngs_course/m2305732_assess/data/reference/hg19.fa \
 --vcf ~/ngs_course/m2305732_assess/results/NGS0001.vcf

#Zip vcf file to save space
bgzip ~/ngs_course/m2305732_assess/results/NGS0001.vcf

tabix -p vcf ~/ngs_course/m2305732_assess/results/NGS0001.vcf.gz

##Quality filter the FreeBayes VCF file
#Install vcflib
conda install vcflib

#Run quality filter by defining filter parameters
vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" \
 ~/ngs_course/m2305732_assess/results/NGS0001.vcf.gz > ~/ngs_course/m2305732_assess/results/NGS0001_filtered.vcf

#Filter VCF to the regions in the bed file
bedtools intersect -header -wa -a ~/ngs_course/m2305732_assess/results/NGS0001_filtered.vcf \
 -b ~/ngs_course/m2305732_assess/data/annotation.bed > ~/ngs_course/m2305732_assess/results/NGS0001_bed_filtered.vcf

#Zip filtered VCF file
bgzip ~/ngs_course/m2305732_assess/results/NGS0001_bed_filtered.vcf

tabix -p vcf ~/ngs_course/m2305732_assess/results/NGS0001_bed_filtered.vcf.gz

##################################################################################################################################

############## STEP 6: ANNOTATION ################################################################################################
##Annotation with Annovar
#Convert VCF to Annovar input format
~/ngs_course/m2305732_assess/annovar/annovar/convert2annovar.pl \
 -format vcf4 ~/ngs_course/m2305732_assess/results/NGS0001_bed_filtered.vcf.gz > \
 ~/ngs_course/m2305732_assess/results/NGS0001_bed_filtered.avinput

#Move into Annovar directory
cd ~/ngs_course/m2305732_assess/annovar/annovar

#Run Annovar table function to output a csv file
./table_annovar.pl ~/ngs_course/m2305732_assess/results/NGS0001_bed_filtered.avinput humandb/ \
 -buildver hg19 -out ~/ngs_course/m2305732_assess/results/NGS0001_bed_filtered -remove \
 -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro -operation g,g,f,f,f -otherinfo -nastring . -csvout

#OUTPUT LOCATION: /home/ubuntu/ngs_course/m2305732_assess/results/NGS0001_bed_filtered.hg19_multianno.csv

##Filter for variants not in dbsnp
./annotate_variation.pl -filter -out ~/ngs_course/m2305732_assess/results/NGS0001 \
 -build hg19 -dbtype snp138 ~/ngs_course/m2305732_assess/results/NGS0001_bed_filtered.avinput humandb/

#OUTPUT LOCATION: /home/ubuntu/ngs_course/m2305732_assess/results/NGS0001.hg19_snp138_filtered
#Above output contains variants not found in dbsnp

#Convert dbsnp filtered file to csv file
./table_annovar.pl ~/ngs_course/m2305732_assess/results/NGS0001.hg19_snp138_filtered humandb/ \
 -buildver hg19 -out ~/ngs_course/m2305732_assess/results/NGS0001.hg19_snp138_filtered -remove \
 -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro -operation g,g,f,f,f -otherinfo -nastring . -csvout

#Transfer csv file to personal device using filezilla, filter in excel for exonic variants.
#Viewable as excel file uploaded to Canvas: "NGS0001.hg19_snp138_filtered.hg19_multianno.csv"

########## Annotation with snpEFF ################################################################################################
##Install snpeff
#Move to results directory
cd ~/ngs_course/m2305732_assess/results

#Download snpEff
wget http://sourceforge.net/projects/snpeff/files/snpEff_v3_6_core.zip

#Unzip snpEff
unzip snpEff_v3_6_core.zip

#Remove zipped snpEff to save space
rm -r snpEff_v3_6_core.zip

#Move into snpeff directory
cd ~/ngs_course/m2305732_assess/results/snpEff

#install hg19 database
java -jar snpEff.jar download hg19

#Unzip hg19 database: hg19 database is unzipped into data directory withing snpEff
unzip snpEff_v3_6_hg19.zip

#Remove zipped hg19 database to save space
rm -r snpEff_v3_6_hg19.zip

#Move back into results directory to access vcf files
cd ~/ngs_course/m2305732_assess/results

##Annotate using Snpeff: "-Xmx4G" defines 4GB memory usage, "-i" & "-o" input & output format is vcf
java -Xmx4G -jar snpEff/snpEff.jar -i vcf -o vcf hg19 NGS0001_bed_filtered.vcf.gz > NGS0001_Snpeff.vcf

#OUTPUT: snpEff_summary.html
#Transfer html to personal device using Filezilla to view

#################################################################################################################################
