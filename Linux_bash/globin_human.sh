##############################################################################
#      RNA-seq analysis of reticulocyte-derived globin gene transcripts      #
#             (HBA and HBB) in human peripheral blood samples                #
#                  --- Linux bioinformatics workflow ---                     #
##############################################################################
# DOI badge: 
# Author: Correia, C.N.
# Version 1.0.0
# Last updated on: 17/02/2016


############################################
# Human: Download raw FASTQ files from ENA #
############################################

# Create and enter working directory:
mkdir $HOME/storage/globin/human_fastq
cd !$

# Get ftp links for FASTQ files from the European Nucleotide Archive (ENA)
# http://www.ebi.ac.uk/ena/data/view/PRJNA232593
# Under the Read Files tab, click on select columns and make sure that the
# only checke option is Fastq files (ftp). Then, on the download option 
# click on the TEXT link and download the PRJNA232593.txt file.
# Transfer it to Stampede using WinSCP.
# Using notepad++, remove the header on PRJNA232593.txt. Then, search ; and
# replace with \n
# Every ftp link should be in one line.

# Download human data set from ENA (PRJNA232593, Shin et al.2014):
nohup wget -i PRJNA232593.txt &

# Check that all files were downloaded:
grep '.fastq.gz' PRJNA232593.txt | wc -l   # 94 lines
grep '.fastq.gz’ saved' nohup.out | wc -l  # 94 lines

##################################################
# Human: FastQC quality check of raw FASTQ files #
##################################################

# Required software is FastQC v0.11.5, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and enter the quality check output directory:
mkdir $HOME/scratch/globin/quality_check/pre-filtering/human
cd !$

# Run FastQC in one file to see if it's working well:
fastqc -o $HOME/scratch/globin/quality_check/pre-filtering/human \
--noextract --nogroup -t 2 \
$HOME/storage/globin/human_fastq/SRR1060753_1.fastq.gz

### Moved this folder to my laptop using WinSCP
### and checked the HTML report. It worked fine.

# Create a bash script to perform FastQC quality check on all fastq.gz files:
for file in `find $HOME/storage/globin/human_fastq/ \
-name *fastq.gz`; do echo "fastqc --noextract --nogroup -t 1 \
-o $HOME/scratch/globin/quality_check/pre-filtering/human $file" \
>> fastqc.sh; done;

# Run script on Stampede:
chmod 755 fastqc.sh
nohup ./fastqc.sh > fastqc.sh.nohup &

# Check if all the files were processed:
ls -l | grep fastqc.zip | wc -l

more fastqc.sh.nohup | grep "Failed to process file" >> failed_fastqc.txt

# Deleted all the HTML files:
rm -r *.html

### Copied all .zip files to my laptop using WinSCP
### and checked the HTML reports. It worked fine.

# Check all output from FastQC:
mkdir $HOME/scratch/globin/quality_check/pre-filtering/human/tmp

for file in `ls *_fastqc.zip`; do unzip \
$file -d $HOME/scratch/globin/quality_check/pre-filtering/human/tmp; \
done

for file in \
`find $HOME/scratch/globin/quality_check/pre-filtering/human/tmp \
-name summary.txt`; do more $file >> reports_pre-filtering.txt; \
done

grep 'Adapter Content' reports_pre-filtering.txt >> adapter_content.txt
wc -l adapter_content.txt

grep FAIL adapter_content.txt | wc -l

for file in \
`find $HOME/scratch/globin/quality_check/pre-filtering/human/tmp \
-name fastqc_data.txt`; do head -n 10 $file >> basic_stats_pre-filtering.txt; \
done

# Remove temporary folder and its files:
rm -r $HOME/scratch/globin/quality_check/pre-filtering/human/tmp

#########################################################################
# Human: Adapter-contamination and quality filtering of raw FASTQ files #
#########################################################################

# Required software is ngsShoRT (version 2.2). More information can be found
# here: http://research.bioinformatics.udel.edu/genomics/ngsShoRT/index.html

# Create a working directory for filtered reads:
mkdir $HOME/scratch/globin/fastq_sequence/human
cd !$

# Copy Illumina adpaters file into working directory:
cp $HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/Illumina_PE_adapters.txt .

# Run ngsShoRT in one pair of reads to check if it's working:
nohup perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 5 -mode trim -min_rl 70 \
-pe1 $HOME/storage/globin/human_fastq/SRR1060753_1.fastq.gz \
-pe2 $HOME/storage/globin/human_fastq/SRR1060753_2.fastq.gz \
-o $HOME/scratch/globin/fastq_sequence/human/SRR1060753 \
-methods 5adpt_lqr -5a_f Illumina_PE_adapters.txt -5a_mp 90 -5a_del 0 \
-5a_ins 0 -5a_fmi 100 -5a_axn kr -lqs 20 -lq_p 25 -gzip &

# Create bash script to perform filtering of each FASTQ file:
for file in `find $HOME/storage/globin/human_fastq \
-name *_1.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/_1.fastq.gz/_2.fastq.gz/'`; \
sample=`basename $file | perl -p -e 's/_1.fastq.gz//'`; \
echo "perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 5 -mode trim -min_rl 70 \
-pe1 $file -pe2 $file2 \
-o $HOME/scratch/globin/fastq_sequence/human/$sample \
-methods 5adpt_lqr -5a_f Illumina_PE_adapters.txt -5a_mp 90 -5a_del 0 \
-5a_ins 0 -5a_fmi 100 -5a_axn kr -lqs 20 -lq_p 25 -gzip" \
>> filtering.sh; \
done;

# Split and run all scripts on Stampede:
split -d -l 16 filtering.sh filtering.sh.
for script in `ls filtering.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Check that all folders were created:
ls -l | grep SRR | wc -l

# Check that all pairs were processed:
for file in `ls filtering.sh.0*.nohup`; \
do grep -o 'Done-MAIN' $file | wc -l; done

# Compress files of removed reads:
for file in `find $HOME/scratch/globin/fastq_sequence/human \
-name extracted*.txt`; \
do echo "gzip -9 $file" >> discarded_compression.sh; \
done

# Run script on Stampede:
chmod 755 discarded_compression.sh
nohup ./discarded_compression.sh &

#######################################################
# Human: FastQC quality check of filtered FASTQ files #
#######################################################

# Required software is FastQC v0.11.5, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and enter the quality check output directory:
mkdir $HOME/scratch/globin/quality_check/post-filtering/human
cd !$

# Run FastQC in one file to see if it's working well:
fastqc -o $HOME/scratch/globin/quality_check/post-filtering/human \
--noextract --nogroup -t 5 \
$HOME/scratch/globin/fastq_sequence/human/SRR1060753/trimmed_SRR1060753_1.fastq.gz

### Moved this folder to my laptop using WinSCP
### and checked the HTML report. It worked fine.

# Create a bash script to perform FastQC quality check on all fastq.gz files:
for file in `find $HOME/scratch/globin/fastq_sequence/human/ \
-name trimmed_*.fastq.gz`; do echo "fastqc --noextract --nogroup -t 1 \
-o $HOME/scratch/globin/quality_check/post-filtering/human $file" \
>> fastqc.sh; done;

# Run script on Stampede:
chmod 755 fastqc.sh
nohup ./fastqc.sh > fastqc.sh.nohup &

# Check if all the files were processed:
ls -l | grep fastqc.zip | wc -l

more fastqc.sh.nohup | grep "Failed to process file" >> failed_fastqc.txt

# Deleted all the HTML files:
rm -r *.html

### Copied all .zip files to my laptop using WinSCP
### and checked the HTML reports. It worked fine.

# Check all output from FastQC:
mkdir $HOME/scratch/globin/quality_check/post-filtering/human/tmp

for file in `ls *_fastqc.zip`; do unzip \
$file -d $HOME/scratch/globin/quality_check/post-filtering/human/tmp; \
done

for file in \
`find $HOME/scratch/globin/quality_check/post-filtering/human/tmp \
-name summary.txt`; do more $file >> reports_post-filtering.txt; \
done

grep 'Adapter Content' reports_post-filtering.txt >> adapter_content.txt
wc -l adapter_content.txt

grep FAIL adapter_content.txt | wc -l

for file in \
`find $HOME/scratch/globin/quality_check/post-filtering/human/tmp \
-name fastqc_data.txt`; do head -n 10 $file >> basic_stats_post-filtering.txt; \
done

# Remove temporary folder and its files:
rm -r $HOME/scratch/globin/quality_check/post-filtering/human/tmp


############################################################
# Human: Alignment of FASTQ files against the Homo sapiens #
#                reference genome with STAR                #
############################################################

### Moved trimmed reads to Rodeo.
### Following steps were conducted in Ubuntu 14.04

# Required software is STAR 2.5.2b, consult manual/tutorial for details:
https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

# Download the latest human reference genome from NCBI RefSeq
# (GCF_000001405.35_GRCh38.p9):
mkdir -p /home/workspace/genomes/homosapiens/hg38_p9/source_file
cd !$

nohup wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.35_GRCh38.p9/GCF_000001405.35_GRCh38.p9_genomic.fna.gz &
gunzip GCF_000001405.35_GRCh38.p9_genomic.fna.gz

# Download annotation file for GCF_000001405.35_GRCh38.p9:
mkdir -p /home/workspace/genomes/homosapiens/hg38_p9/annotation_file
cd !$

nohup wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.35_GRCh38.p9/GCF_000001405.35_GRCh38.p9_genomic.gff.gz &
gunzip GCF_000001405.35_GRCh38.p9_genomic.gff.gz

# Generate the genome index using annotations:
mkdir /home/workspace/genomes/homosapiens/hg38_p9/STAR-2.5.2b_index_74bp
cd !$

nohup STAR --runThreadN 40 --runMode genomeGenerate \
--genomeDir /home/workspace/genomes/homosapiens/hg38_p9/STAR-2.5.2b_index_74bp \
--genomeFastaFiles \
/home/workspace/genomes/homosapiens/hg38_p9/source_file/GCF_000001405.35_GRCh38.p9_genomic.fna \
--sjdbGTFfile /home/workspace/genomes/homosapiens/hg38_p9/annotation_file/GCF_000001405.35_GRCh38.p9_genomic.gff \
--sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 74 \
--outFileNamePrefix \
/home/workspace/genomes/homosapiens/hg38_p9/STAR-2.5.2b_index_74bp/h38.p9 &

# Create and enter alignment working directory:
mkdir -p /home/workspace/ccorreia/globin/STAR-2.5.2b_alignment/human
cd !$

# Mapping reads from one FASTQ file to the indexed genome,
# to check if it works well:
STAR --runMode alignReads --runThreadN 20 --genomeLoad LoadAndRemove \
--genomeDir /home/workspace/genomes/homosapiens/hg38_p9/STAR-2.5.2b_index_74bp \
--readFilesIn \
/home/workspace/ccorreia/globin/fastq_sequence/human/SRR1060753/trimmed_SRR1060753_1.fastq.gz \
/home/workspace/ccorreia/globin/fastq_sequence/human/SRR1060753/trimmed_SRR1060753_2.fastq.gz \
--readFilesCommand gunzip -c --outFilterMultimapNmax 10 \
--outFilterMismatchNmax 10 --outFileNamePrefix ./SRR1060753_ \
--outSAMtype BAM Unsorted --outReadsUnmapped Fastx

# Create a bash script to perform alignment of paired FASTQ files:
for file in `find /home/workspace/ccorreia/globin/fastq_sequence/human \
-name *_1.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/\_1\.fastq\.gz/\_2\.fastq\.gz/'`; \
sample=`basename $file | perl -p -e 's/\_1\.fastq\.gz//'`; \
echo "mkdir /home/workspace/ccorreia/globin/STAR-2.5.2b_alignment/human/$sample; \
cd /home/workspace/ccorreia/globin/STAR-2.5.2b_alignment/human/$sample; \
STAR --runMode alignReads --runThreadN 20 --genomeLoad LoadAndRemove \
--genomeDir /home/workspace/genomes/homosapiens/hg38_p9/STAR-2.5.2b_index_74bp \
--readFilesIn $file $file2 --readFilesCommand gunzip -c \
--outFilterMultimapNmax 10 --outFilterMismatchNmax 10 \
--outFileNamePrefix ./${sample}_ --outSAMtype BAM Unsorted \
--outSAMattrIHstart 0 --outSAMattributes Standard --outReadsUnmapped Fastx" \
>> alignment.sh; \
done

# Split and run script on Rodeo:
split -d -l 16 alignment.sh alignment.sh.
for script in `ls alignment.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Check nohup.out file to see how many jobs finished successfully:
grep -c 'finished successfully' alignment.sh.00.nohup
grep -c 'finished successfully' alignment.sh.01.nohup
grep -c 'finished successfully' alignment.sh.02.nohup

# Merge all STAR log.final.out files into a single file:
for file in `find /home/workspace/ccorreia/globin/STAR-2.5.2b_alignment/human \
-name *Log.final.out`; \
do perl /home/workspace/ccorreia/scripts/star_report_opener.pl -report $file; \
done

#############################################
# FastQC quality check of aligned BAM files #
#############################################

# Required software is FastQC v0.11.5, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and go to working directory:
mkdir -p /home/workspace/ccorreia/globin/quality_check/post_alignment/human
cd !$

# Create a bash script to perform FastQC quality check on aligned BAM files:
for file in `find /home/workspace/ccorreia/globin/STAR-2.5.2b_alignment/human \
-name *.bam`; do echo "fastqc-0.11.5 --noextract --nogroup -t 10 \
-o /home/workspace/ccorreia/globin/quality_check/post_alignment/human $file" >> \
fastqc_aligned.sh; \
done

# Run script on Rodeo:
chmod 755 fastqc_aligned.sh
nohup ./fastqc_aligned.sh > fastqc_aligned.sh.nohup &

# Delete all the HTML files:
rm -r *.html

# Check all output from FastQC:
mkdir /home/workspace/ccorreia/globin/quality_check/post_alignment/human/tmp

for file in `ls *_fastqc.zip`; do unzip \
$file -d /home/workspace/ccorreia/globin/quality_check/post_alignment/human/tmp; \
done

for file in \
`find /home/workspace/ccorreia/globin/quality_check/post_alignment/human/tmp \
-name summary.txt`; do more $file >> reports_post-alignment.txt; \
done

for file in \
`find /home/workspace/ccorreia/globin/quality_check/post_alignment/human/tmp \
-name fastqc_data.txt`; do head -n 10 $file >> basic_stats_post_alignment.txt; \
done

# Check if all files were processed:
grep -c '##FastQC' basic_stats_post_alignment.txt
grep -c 'Basic Statistics' reports_post-alignment.txt
grep -c 'Analysis complete' fastqc_aligned.sh.nohup

# Remove temporary folder:
rm -r tmp/

###################################################################
# Summarisation of gene counts with featureCounts for sense genes #
###################################################################

# Required package is featureCounts, which is part of Subread 1.5.1 software,
# consult manual for details:
# http://bioinf.wehi.edu.au/subread-package/SubreadUsersGuide.pdf

# Create working directories:
mkdir -p /home/workspace/ccorreia/globin/Count_summarisation/sense/human
cd !$

# Run featureCounts with one sample to check if it is working fine:
featureCounts -a \
/home/workspace/genomes/homosapiens/hg38_p9/annotation_file/GCF_000001405.35_GRCh38.p9_genomic.gff \
-B -p -C -R -s 1 -T 20 -t gene -g Dbxref -o ./counts.txt \
/home/workspace/ccorreia/globin/STAR-2.5.2b_alignment/human/trimmed_SRR1060753/trimmed_SRR1060753_Aligned.out.bam

# Create a bash script to run featureCounts on BAM file containing multihits and
# uniquely mapped reads using the stranded parameter:
for file in `find /home/workspace/ccorreia/globin/STAR-2.5.2b_alignment/human \
-name *_Aligned.out.bam`; \
do sample=`basename $file | perl -p -e 's/_Aligned.out.bam//'`; \
echo "mkdir /home/workspace/ccorreia/globin/Count_summarisation/sense/human/$sample; \
cd /home/workspace/ccorreia/globin/Count_summarisation/sense/human/$sample; \
featureCounts -a \
/home/workspace/genomes/homosapiens/hg38_p9/annotation_file/GCF_000001405.35_GRCh38.p9_genomic.gff \
-B -p -C -R -s 1 -T 20 -t gene -g Dbxref \
-o ${sample}_sense-counts.txt $file" >> sense_count.sh; \
done

# Split and run all scripts on Rodeo:
split -d -l 16 sense_count.sh sense_count.sh.
for script in `ls sense_count.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Check if all files were processed:
grep -c 'Read assignment finished.' sense_count.sh.00.nohup
grep -c 'Read assignment finished.' sense_count.sh.01.nohup
grep -c 'Read assignment finished.' sense_count.sh.02.nohup

# Create bash script to merge stats info from .featureCounts from all samples
# into a single file:
for file in `find /home/workspace/ccorreia/globin/Count_summarisation/sense/human \
-name *.featureCounts`; do echo echo \
"\`basename $file\` \`cut $file -f2 | sort | uniq -c | perl -p -e 's/\n/ /'\` >> \
annotation_summary_sense.txt" >> annotation_summary_sense.sh
done

# Split and run all scripts on Rodeo:
split -d -l 10 annotation_summary_sense.sh annotation_summary_sense.sh.
for script in `ls annotation_summary_sense.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Check that all files were processed:
grep -c '.featureCounts' annotation_summary_sense.txt

# Copy all *sense-counts.txt files to temporary folder:
mkdir /home/workspace/ccorreia/globin/Count_summarisation/sense/human/tmp

for file in `find /home/workspace/ccorreia/globin/Count_summarisation/sense/human \
-name *sense-counts.txt`; do cp $file \
-t /home/workspace/ccorreia/globin/Count_summarisation/sense/human/tmp; \
done

# Transfer all files from tmp to laptop, using WinSCP, then remove tmp folder:
rm -r tmp



















