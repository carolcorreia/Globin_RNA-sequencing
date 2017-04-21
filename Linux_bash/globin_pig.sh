##############################################################################
#      RNA-seq analysis of reticulocyte-derived globin gene transcripts      #
#         (HBA and HBB) in pig (Sus scrofa) peripheral blood samples         #
#                  --- Linux bioinformatics workflow ---                     #
##############################################################################

# DOI badge: 
# Author: Correia, C.N.
# Version 1.0.0
# Last updated on: 18/04/2017


#################################
# Pig: Download raw FASTQ files #
#################################

# Create and enter working directory:
mkdir $HOME/storage/globin/pig_fastq
cd !$

# Download pig data set as per authors' instructions
# (personal communication, not shown here).
# Choi, I, Bao, H, Kommadath, A, Hosseini, A, Sun, X, Meng, Y, Stothard, P,
# Plastow, GS, Tuggle, CK, Reecy, JM, Fritz-Waters, E, Abrams, SM, Lunney, JK,
# and Guan le, L (2014).
# Increasing gene discovery and coverage using RNA-seq of globin RNA reduced
# porcine blood samples. BMC genomics 15, 954. doi: 10.1186/1471-2164-15-954.

# Check that all files were downloaded:
ls -l | grep fastq.gz | wc -l              # Result: 80
grep '.fastq.gz’ saved' nohup.out | wc -l  # Result: 80

################################################
# Pig: FastQC quality check of raw FASTQ files #
################################################

# Required software is FastQC v0.11.5, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and enter the quality check output directory:
mkdir -p $HOME/scratch/globin/quality_check/pre-filtering/pig
cd !$

# Run FastQC in one file to see if it's working well:
fastqc -o $HOME/scratch/globin/quality_check/pre-filtering/pig \
--noextract --nogroup -t 2 \
$HOME/storage/globin/pig_fastq/HI.0751.004.Index_12.GCswine-5037-28DPI-WB-7413C-mRNA_R1.fastq.gz

### Moved this folder to my laptop using WinSCP
### and checked the HTML report. It worked fine.

# Create a bash script to perform FastQC quality check on all fastq.gz files:
for file in `find $HOME/storage/globin/pig_fastq/ \
-name *fastq.gz`; do echo "fastqc --noextract --nogroup -t 1 \
-o $HOME/scratch/globin/quality_check/pre-filtering/pig $file" \
>> fastqc.sh; done;

# Split and run all scripts on Stampede:
split -d -l 40 fastqc.sh fastqc.sh.
for script in `ls fastqc.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Check if all the files were processed:
ls -l | grep fastqc.zip | wc -l

for file in `ls fastqc.sh.0*.nohup`; \
do more $file | grep "Failed to process file" >> failed_fastqc.txt
done

# Deleted all the HTML files:
rm -r *.html

### Copied all .zip files to my laptop using WinSCP
### and checked the HTML reports. It worked fine.

# Check all output from FastQC:
mkdir $HOME/scratch/globin/quality_check/pre-filtering/pig/tmp

for file in `ls *_fastqc.zip`; do unzip \
$file -d $HOME/scratch/globin/quality_check/pre-filtering/pig/tmp; \
done;

for file in \
`find $HOME/scratch/globin/quality_check/pre-filtering/pig/tmp \
-name summary.txt`; do more $file >> reports_pre-filtering.txt; done

grep 'Adapter Content' reports_pre-filtering.txt >> adapter_content.txt
wc -l adapter_content.txt

for file in \
`find $HOME/scratch/globin/quality_check/pre-filtering/pig/tmp \
-name fastqc_data.txt`; do head -n 10 $file >> basic_stats_pre-filtering.txt; \
done

# Remove temporary folder and its files:
rm -r $HOME/scratch/globin/quality_check/pre-filtering/pig/tmp

####################################
# Pig: Trimming of raw FASTQ files #
####################################

# Required software is ngsShoRT (version 2.2). More information can be found
# here: http://research.bioinformatics.udel.edu/genomics/ngsShoRT/index.html

# Create a working directory for filtered reads:
mkdir $HOME/scratch/globin/fastq_sequence/pig
cd !$

# Copy Illumina adpaters file into working directory:
cp $HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/Illumina_PE_adapters.txt .

# Run ngsShoRT in one pair of reads to check if it's working:
nohup perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 7 -mode trim -min_rl 80 \
-pe1 $HOME/storage/globin/pig_fastq/HI.0751.004.Index_12.GCswine-5037-28DPI-WB-7413C-mRNA_R1.fastq.gz \
-pe2 $HOME/storage/globin/pig_fastq/HI.0751.004.Index_12.GCswine-5037-28DPI-WB-7413C-mRNA_R2.fastq.gz \
-o $HOME/scratch/globin/fastq_sequence/pig \
-methods 5adpt_lqr_3end -5a_f Illumina_PE_adapters.txt -5a_mp 90 -5a_del 0 \
-5a_ins 0 -5a_fmi 100 -5a_axn kr -lqs 20 -lq_p 25 -n3 10 -gzip &

# Create bash scripts to perform trimming of reads (10bp at 3' end), while 
# keeping the sequencing lane information:
for file in `find $HOME/storage/globin/pig_fastq \
-name *_R1.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/_R1.fastq.gz/_R2.fastq.gz/'`; \
sample=`basename $file | perl -p -e 's/-mRNA_R1.fastq.gz//'`; \
echo "perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 5 -mode trim -min_rl 80 \
-pe1 $file -pe2 $file2 \
-o $HOME/scratch/globin/fastq_sequence/pig/$sample \
-methods 5adpt_lqr_3end -5a_f Illumina_PE_adapters.txt -5a_mp 90 -5a_del 0 \
-5a_ins 0 -5a_fmi 100 -5a_axn kr -lqs 20 -lq_p 25 -n3 10 -gzip" \
>> filtering.sh; \
done;

# Split and run all scripts on Stampede:
split -d -l 20 filtering.sh filtering.sh.
for script in `ls filtering.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Check that all folders were created:
ls -l | grep HI | wc -l

# Check that all the pairs were processed:
for file in `ls filtering.sh.0*.nohup`; \
do grep -o 'Done-MAIN' $file | wc -l; done

# Compress files of removed reads:
for file in `find $HOME/scratch/globin/fastq_sequence/pig \
-name extracted*.txt`; \
do echo "gzip -9 $file" >> discarded_compression.sh; \
done

# Run script on Stampede:
chmod 755 discarded_compression.sh
nohup ./discarded_compression.sh &

####################################################
# Pig: FastQC quality check of trimmed FASTQ files #
####################################################

# Required software is FastQC v0.11.5, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and enter the quality check output directory:
mkdir -p $HOME/scratch/globin/quality_check/post-filtering/pig
cd !$

# Run FastQC in one file to see if it's working well:
fastqc -o $HOME/scratch/globin/quality_check/post-filtering/pig \
--noextract --nogroup -t 5 \
$HOME/scratch/globin/fastq_sequence/pig/HI.0751.004.Index_12.GCswine-5037-28DPI-WB-7413C/trimmed_HI.0751.004.Index_12.GCswine-5037-28DPI-WB-7413C-mRNA_R1.fastq.gz

### Moved this folder to my laptop using WinSCP
### and checked the HTML report. It worked fine.

# Create a bash script to perform FastQC quality check on all fastq.gz files:
for file in `find $HOME/scratch/globin/fastq_sequence/pig \
-name trimmed_*fastq.gz`; do echo "fastqc --noextract --nogroup -t 1 \
-o $HOME/scratch/globin/quality_check/post-filtering/pig $file" \
>> fastqc.sh; \
done;

# Split and run all scripts on Stampede:
split -d -l 40 fastqc.sh fastqc.sh.
for script in `ls fastqc.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Check if all the files were processed:
ls -l | grep fastqc.zip | wc -l

for file in `ls fastqc.sh.0*.nohup`; \
do more $file | grep "Failed to process file" >> failed_fastqc.txt
done

# Deleted all the HTML files:
rm -r *.html

### Copied all .zip files to my laptop using WinSCP
### and checked the HTML reports. It worked fine.

# Check all output from FastQC:
mkdir $HOME/scratch/globin/quality_check/post-filtering/pig/tmp

for file in `ls *_fastqc.zip`; do unzip \
$file -d $HOME/scratch/globin/quality_check/post-filtering/pig/tmp; \
done;

for file in \
`find $HOME/scratch/globin/quality_check/post-filtering/pig/tmp \
-name summary.txt`; do more $file >> reports_post-filtering.txt; done

grep 'Adapter Content' reports_post-filtering.txt >> adapter_content.txt
wc -l adapter_content.txt

grep FAIL adapter_content.txt | wc -l

for file in \
`find $HOME/scratch/globin/quality_check/post-filtering/pig/tmp \
-name fastqc_data.txt`; do head -n 10 $file >> basic_stats_post-filtering.txt; \
done

# Remove temporary folder and its files:
rm -r $HOME/scratch/globin/quality_check/post-filtering/pig/tmp

########################################################
# Pig: Alignment of FASTQ files against the Sus scrofa #
#            reference genome with STAR                #
########################################################

### Moved trimmed reads to Rodeo.
### Following steps were conducted in Ubuntu 14.04.

# Required software is STAR 2.5.2b, consult manual/tutorial for details:
https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

# Download Sus scrofa reference genome, version 10.2 from NCBI RefSeq:
mkdir -p /home/workspace/genomes/susscrofa/NCBI_10.2/source_file
cd !$

nohup wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/025/GCF_000003025.5_Sscrofa10.2/GCF_000003025.5_Sscrofa10.2_genomic.fna.gz &
gunzip GCF_000003025.5_Sscrofa10.2_genomic.fna.gz

# Download annotation file for Sscrofa10.2 NCBI RefSeq:
mkdir /home/workspace/genomes/susscrofa/NCBI_10.2/annotation_file
cd !$

nohup wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/025/GCF_000003025.5_Sscrofa10.2/GCF_000003025.5_Sscrofa10.2_genomic.gff.gz &
gunzip GCF_000003025.5_Sscrofa10.2_genomic.gff.gz

# Generate the genome index using annotations:
mkdir /home/workspace/genomes/susscrofa/NCBI_10.2/STAR-2.5.2b_index_89bp
cd !$

nohup STAR --runThreadN 40 --runMode genomeGenerate \
--genomeDir /home/workspace/genomes/susscrofa/NCBI_10.2/STAR-2.5.2b_index_89bp \
--genomeFastaFiles \
/home/workspace/genomes/susscrofa/NCBI_10.2/source_file/GCF_000003025.5_Sscrofa10.2_genomic.fna \
--sjdbGTFfile /home/workspace/genomes/susscrofa/NCBI_10.2/annotation_file/GCF_000003025.5_Sscrofa10.2_genomic.gff \
--sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 89 \
--outFileNamePrefix \
/home/workspace/genomes/susscrofa/NCBI_10.2/STAR-2.5.2b_index_89bp/h38.p9 &

# Create and enter alignment working directory:
mkdir -p /home/workspace/ccorreia/globin/STAR-2.5.2b_alignment/pig
cd !$

# Mapping reads from one FASTQ file to the indexed genome,
# to check if it works well:
STAR --runMode alignReads --runThreadN 20 --genomeLoad LoadAndRemove \
--genomeDir /home/workspace/genomes/susscrofa/NCBI_10.2/STAR-2.5.2b_index_89bp \
--readFilesIn \
/home/workspace/ccorreia/globin/fastq_sequence/pig/HI.0751.004.Index_12.GCswine-5037-28DPI-WB-7413C/trimmed_HI.0751.004.Index_12.GCswine-5037-28DPI-WB-7413C-mRNA_R1.fastq.gz \
/home/workspace/ccorreia/globin/fastq_sequence/pig/HI.0751.004.Index_12.GCswine-5037-28DPI-WB-7413C/trimmed_HI.0751.004.Index_12.GCswine-5037-28DPI-WB-7413C-mRNA_R2.fastq.gz \
--readFilesCommand gunzip -c --outFilterMultimapNmax 10 \
--outFilterMismatchNmax 10 --outFileNamePrefix ./HI.0751.004.Index_12.GCswine-5037-28DPI-WB-7413C_ \
--outSAMtype BAM Unsorted --outReadsUnmapped Fastx

# Create a bash script to perform alignment of paired FASTQ files sequenced
# over different lanes.
# Lanes 002 or 003:
for file in `find /home/workspace/ccorreia/globin/fastq_sequence/pig \
-name trimmed_*.00[23].*_R1.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/\_R1\.fastq\.gz/\_R2\.fastq\.gz/'`; \
sample=`basename $file | perl -p -e 's/trimmed_.*(7\d\d\d.+)\-mRNA_R1\.fastq\.gz/$1/'`; \
echo "mkdir /home/workspace/ccorreia/globin/STAR-2.5.2b_alignment/pig/$sample; \
cd /home/workspace/ccorreia/globin/STAR-2.5.2b_alignment/pig/$sample; \
STAR --runMode alignReads --runThreadN 20 --genomeLoad LoadAndRemove \
--genomeDir /home/workspace/genomes/susscrofa/NCBI_10.2/STAR-2.5.2b_index_89bp \
--readFilesIn $file $file2 --readFilesCommand gunzip -c \
--outFilterMultimapNmax 10 --outFilterMismatchNmax 10 \
--outFileNamePrefix ./${sample}_ --outSAMtype BAM Unsorted \
--outSAMattrIHstart 0 --outSAMattributes Standard --outReadsUnmapped Fastx" \
>> alignment.sh; \
done

# Lanes 004 and 005:
for file in `find /home/workspace/ccorreia/globin/fastq_sequence/pig \
-name trimmed_*.004.*_R1.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/(\.004\.)/\.005\./g'`; \
read1=`echo $file | perl -p -e 's/\_R1\.fastq\.gz/\_R2\.fastq\.gz/'`; \
read2=`echo $file2 | perl -p -e 's/\_R1\.fastq\.gz/\_R2\.fastq\.gz/'`; \
sample=`basename $file | perl -p -e 's/trimmed_.*(7\d\d\d.+)\-mRNA_R1\.fastq\.gz/$1/'`; \
echo "mkdir /home/workspace/ccorreia/globin/STAR-2.5.2b_alignment/pig/$sample; \
cd /home/workspace/ccorreia/globin/STAR-2.5.2b_alignment/pig/$sample; \
STAR --runMode alignReads --runThreadN 20 --genomeLoad LoadAndRemove \
--genomeDir /home/workspace/genomes/susscrofa/NCBI_10.2/STAR-2.5.2b_index_89bp \
--readFilesIn $file,$file2 $read1,$read2 --readFilesCommand gunzip -c \
--outFilterMultimapNmax 10 --outFilterMismatchNmax 10 \
--outFileNamePrefix ./${sample}_ --outSAMtype BAM Unsorted \
--outSAMattrIHstart 0 --outSAMattributes Standard --outReadsUnmapped Fastx" \
>> alignment.sh; \
done

# Split and run all scripts on Rodeo:
split -d -l 12 alignment.sh alignment.sh.
for script in `ls alignment.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Check nohup.out file to see how many jobs finished successfully:
grep -c 'finished successfully' alignment.sh.00.nohup
grep -c 'finished successfully' alignment.sh.01.nohup

# Merge all STAR log.final.out files into a single file:
for file in `find /home/workspace/ccorreia/globin/STAR-2.5.2b_alignment/pig \
-name *Log.final.out`; \
do perl /home/workspace/ccorreia/scripts/star_report_opener.pl -report $file; \
done

#############################################
# FastQC quality check of aligned BAM files #
#############################################

# Required software is FastQC v0.11.5, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and go to working directory:
mkdir -p /home/workspace/ccorreia/globin/quality_check/post_alignment/pig
cd !$

# Create a bash script to perform FastQC quality check on aligned BAM files:
for file in `find /home/workspace/ccorreia/globin/STAR-2.5.2b_alignment/pig \
-name *.bam`; do echo "fastqc --noextract --nogroup -t 10 \
-o /home/workspace/ccorreia/globin/quality_check/post_alignment/pig $file" >> \
fastqc_aligned.sh; \
done

# Split and run all scripts on Rodeo:
split -d -l 6 fastqc_aligned.sh fastqc_aligned.sh.
for script in `ls fastqc_aligned.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Delete all the HTML files:
rm -r *.html

# Check all output from FastQC:
mkdir /home/workspace/ccorreia/globin/quality_check/post_alignment/pig/tmp

for file in `ls *_fastqc.zip`; do unzip \
$file -d /home/workspace/ccorreia/globin/quality_check/post_alignment/pig/tmp; \
done

for file in \
`find /home/workspace/ccorreia/globin/quality_check/post_alignment/pig/tmp \
-name summary.txt`; do more $file >> reports_post-alignment.txt; \
done

for file in \
`find /home/workspace/ccorreia/globin/quality_check/post_alignment/pig/tmp \
-name fastqc_data.txt`; do head -n 10 $file >> basic_stats_post_alignment.txt; \
done

# Check if all files were processed:
grep -c '##FastQC' basic_stats_post_alignment.txt
grep -c 'Basic Statistics' reports_post-alignment.txt
grep -c 'Analysis complete' fastqc_aligned.sh.00.nohup
grep -c 'Analysis complete' fastqc_aligned.sh.01.nohup
grep -c 'Analysis complete' fastqc_aligned.sh.02.nohup
grep -c 'Analysis complete' fastqc_aligned.sh.03.nohup

# Transfer FASTQC files to laptop:
scp -r remote_server:/home/workspace/ccorreia/globin/quality_check/post_alignment/pig/tmp .

# Remove tmp folder from Rodeo:
rm -r tmp

###################################################################
# Summarisation of gene counts with featureCounts for sense genes #
###################################################################

# Required package is featureCounts, which is part of Subread 1.5.1 software,
# consult manual for details:
# http://bioinf.wehi.edu.au/subread-package/SubreadUsersGuide.pdf

# Create working directories:
mkdir -p /home/workspace/ccorreia/globin/Count_summarisation/sense/pig
cd !$

# Run featureCounts with one sample to check if it is working fine:
featureCounts -a \
/home/workspace/genomes/susscrofa/NCBI_10.2/annotation_file/GCF_000003025.5_Sscrofa10.2_genomic.gff \
-B -p -C -R -s 0 -T 20 -t gene -g Dbxref -o ./unstranded_counts.txt \
/home/workspace/ccorreia/globin/STAR-2.5.2b_alignment/pig/7413C/7413C_Aligned.out.bam

# Create a bash script to run featureCounts on BAM file containing multihits and
# uniquely mapped reads using the stranded parameter:
for file in `find /home/workspace/ccorreia/globin/STAR-2.5.2b_alignment/pig \
-name *_Aligned.out.bam`; \
do sample=`basename $file | perl -p -e 's/_Aligned.out.bam//'`; \
echo "mkdir /home/workspace/ccorreia/globin/Count_summarisation/sense/pig/$sample; \
cd /home/workspace/ccorreia/globin/Count_summarisation/sense/pig/$sample; \
featureCounts -a \
/home/workspace/genomes/susscrofa/NCBI_10.2/annotation_file/GCF_000003025.5_Sscrofa10.2_genomic.gff \
-B -p -C -R -s 0 -T 20 -t gene -g Dbxref \
-o ${sample}_unstranded-counts.txt $file" >> unstranded_count.sh; \
done

# Split and run all scripts on Rodeo:
split -d -l 12 unstranded_count.sh unstranded_count.sh.
for script in `ls unstranded_count.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Check if all files were processed:
grep -c 'Read assignment finished.' unstranded_count.sh.00.nohup
grep -c 'Read assignment finished.' unstranded_count.sh.01.nohup

# Create bash script to merge stats info from .featureCounts from all samples
# into a single file:
for file in `find /home/workspace/ccorreia/globin/Count_summarisation/sense/pig \
-name *.featureCounts`; do echo echo \
"\`basename $file\` \`cut $file -f2 | sort | uniq -c | perl -p -e 's/\n/ /'\` >> \
annotation_summary.txt" >> annotation_summary.sh
done

# Split and run all scripts on Rodeo:
split -d -l 6 annotation_summary.sh annotation_summary.sh.
for script in `ls annotation_summary.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Check that all files were processed:
grep -c '.featureCounts' annotation_summary.txt

# Copy all *sense-counts.txt files to temporary folder:
mkdir /home/workspace/ccorreia/globin/Count_summarisation/sense/pig/tmp

for file in `find /home/workspace/ccorreia/globin/Count_summarisation/sense/pig \
-name *unstranded-counts.txt`; do cp $file \
-t /home/workspace/ccorreia/globin/Count_summarisation/sense/pig/tmp; \
done

# Transfer all files from tmp to laptop, then remove tmp folder:
scp -r remote_server@rodeo.ucd.ie:/home/workspace/ccorreia/globin/Count_summarisation/sense/pig/tmp .

# Remove tmp folder from Rodeo:
rm -r tmp

# Following steps will be performed in R.















