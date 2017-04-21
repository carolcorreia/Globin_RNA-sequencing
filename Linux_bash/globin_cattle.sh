##############################################################################
#      RNA-seq analysis of reticulocyte-derived globin gene transcripts      #
#             (HBA and HBB) in cattle peripheral blood samples               #
#                  --- Linux bioinformatics workflow ---                     #
##############################################################################

# DOI badge: 
# Author: Correia, C.N.
# Version 1.0.0
# Last updated on: 20/04/2017

mkdir /home/ccorreia/scratch/PPDbRNAseqTimeCourse/fastq_sequence/tmp005

for file in `ls /home/workspace/ccorreia/globin/fastq_sequence/cattle/tmp005/*.fastq.gz`; \
do cp $file \
-t /home/workspace/ccorreia/globin/fastq_sequence/cattle; \
done

##############################################################################
# Alignment of FASTQ files against the Bos taurus reference genome with STAR #
##############################################################################

# Required software is STAR 2.5.1b, consult manual/tutorial for details:
https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

# Download Bos taurus reference genome, version UMD3.1.1 from NCBI:
mkdir -p /home/workspace/genomes/bostaurus/UMD3.1.1_NCBI/source_file
cd !$

nohup wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/055/GCF_000003055.6_Bos_taurus_UMD_3.1.1/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.fna.gz &
gunzip GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.fna.gz

# Download annotation file for UMD3.1.1 NCBI Bos taurus Annotation Release 105:
mkdir /home/workspace/genomes/bostaurus/UMD3.1.1_NCBI/annotation_file
cd !$

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/055/GCF_000003055.6_Bos_taurus_UMD_3.1.1/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.gff.gz
gunzip GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.gff.gz

# Generate genome indexes files using annotations:
mkdir /home/workspace/genomes/bostaurus/UMD3.1.1_NCBI/STAR-2.5.2b_index_99
cd !$

nohup STAR --runThreadN 30 --runMode genomeGenerate \
--genomeDir /home/workspace/genomes/bostaurus/UMD3.1.1_NCBI/STAR-2.5.2b_index_99 \
--genomeFastaFiles \
/home/workspace/genomes/bostaurus/UMD3.1.1_NCBI/source_file/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.fna \
--sjdbGTFfile /home/workspace/genomes/bostaurus/UMD3.1.1_NCBI/annotation_file/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.gff \
--sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 99 \
--outFileNamePrefix \
/home/workspace/genomes/bostaurus/UMD3.1.1_NCBI/STAR-2.5.2b_index_99/Bt-UMD3.1.1 &

# Create and enter alignment working directory:
mkdir /home/workspace/ccorreia/globin/STAR-2.5.2b_alignment/cattle
cd !$

# Mapping reads from one FASTQ file to the indexed genome,
# to check if it works well:
nohup STAR --runMode alignReads --runThreadN 20 --genomeLoad LoadAndRemove \
--genomeDir /home/workspace/genomes/bostaurus/UMD3.1.1_NCBI/STAR-2.5.2b_index_99 \
--readFilesIn \
/home/workspace/ccorreia/globin/fastq_sequence/cattle/trimmed_A6511_W-1_F_R1_001.fastq.gz,\
/home/workspace/ccorreia/globin/fastq_sequence/cattle/trimmed_A6511_W-1_F_R1_002.fastq.gz,\
/home/workspace/ccorreia/globin/fastq_sequence/cattle/trimmed_A6511_W-1_F_R1_003.fastq.gz,\
/home/workspace/ccorreia/globin/fastq_sequence/cattle/trimmed_A6511_W-1_F_R1_004.fastq.gz,\
/home/workspace/ccorreia/globin/fastq_sequence/cattle/trimmed_A6511_W-1_F_R1_005.fastq.gz \
/home/workspace/ccorreia/globin/fastq_sequence/cattle/trimmed_A6511_W-1_F_R2_001.fastq.gz,\
/home/workspace/ccorreia/globin/fastq_sequence/cattle/trimmed_A6511_W-1_F_R2_002.fastq.gz,\
/home/workspace/ccorreia/globin/fastq_sequence/cattle/trimmed_A6511_W-1_F_R2_003.fastq.gz,\
/home/workspace/ccorreia/globin/fastq_sequence/cattle/trimmed_A6511_W-1_F_R2_004.fastq.gz,\
/home/workspace/ccorreia/globin/fastq_sequence/cattle/trimmed_A6511_W-1_F_R2_005.fastq.gz \
--readFilesCommand gunzip -c --outFilterMultimapNmax 10 \
--outFilterMismatchNmax 10 --outFileNamePrefix ./A6511_W10_F_ \
--outSAMtype BAM Unsorted --outReadsUnmapped Fastx &

# Create a bash script to perform alignment of paired FASTQ files:
for file in `find /home/workspace/ccorreia/globin/fastq_sequence/cattle \
-name *_R1_001.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/(_00.)/_002/g'`; \
file3=`echo $file | perl -p -e 's/(_00.)/_003/g'`; \
file4=`echo $file | perl -p -e 's/(_00.)/_004/g'`; \
file5=`echo $file | perl -p -e 's/(_00.)/_005/g'`; \
read1=`echo $file | perl -p -e 's/(R1_00.)/R2_001/'`; \
read2=`echo $file2 | perl -p -e 's/(R1_00.)/R2_002/'`; \
read3=`echo $file3 | perl -p -e 's/(R1_00.)/R2_003/'`; \
read4=`echo $file4 | perl -p -e 's/(R1_00.)/R2_004/'`; \
read5=`echo $file5 | perl -p -e 's/(R1_00.)/R2_005/'`; \
sample=`basename $file | perl -p -e 's/trimmed_(A\d\d\d\d_W\-1_F).*\.fastq\.gz/$1/'`; \
echo "mkdir /home/workspace/ccorreia/globin/STAR-2.5.2b_alignment/cattle/$sample; \
cd /home/workspace/ccorreia/globin/STAR-2.5.2b_alignment/cattle/$sample; \
STAR --runMode alignReads --runThreadN 1 --genomeLoad LoadAndRemove \
--genomeDir /home/workspace/genomes/bostaurus/UMD3.1.1_NCBI/STAR-2.5.2b_index_99 \
--readFilesIn $file,$file2,$file3,$file4,$file5 \
$read1,$read2,$read3,$read4,$read5 --readFilesCommand gunzip -c \
--outFilterMultimapNmax 10 --outFilterMismatchNmax 10 \
--outFileNamePrefix ./${sample}_ --outSAMtype BAM Unsorted \
--outSAMattrIHstart 0 --outSAMattributes Standard --outReadsUnmapped Fastx" \
>> alignment.sh; \
done

# Split and run all scripts on Rodeo:
split -d -l 5 alignment.sh alignment.sh.
for script in `ls alignment.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Check nohup.out file to see how many jobs finished successfully:
grep -c 'finished successfully' alignment.sh.00.nohup
grep -c 'finished successfully' alignment.sh.01.nohup

# Merge all STAR log.final.out files into a single file:
for file in `find /home/workspace/ccorreia/globin/STAR-2.5.2b_alignment/cattle \
-name *Log.final.out`; \
do perl /home/workspace/ccorreia/scripts/star_report_opener.pl -report $file; \
done

#############################################
# FastQC quality check of aligned BAM files #
#############################################

# Required software is FastQC v0.11.5, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and go to working directory:
mkdir /home/workspace/ccorreia/globin/quality_check/post_alignment/cattle
cd !$

# Create a bash script to perform FastQC quality check on aligned SAM files:
for file in `find /home/workspace/ccorreia/globin/STAR-2.5.2b_alignment/cattle \
-name *.bam`; do echo "fastqc --noextract --nogroup -t 10 \
-o /home/workspace/ccorreia/globin/quality_check/post_alignment/cattle \
$file" >> fastqc_aligned.sh; \
done

# Split and run all scripts on Rodeo:
split -d -l 5 fastqc_aligned.sh fastqc_aligned.sh.
for script in `ls fastqc_aligned.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Delete all the HTML files:
rm -r *.html

# Check all output from FastQC:
mkdir /home/workspace/ccorreia/globin/quality_check/post_alignment/cattle/tmp

for file in `ls *_fastqc.zip`; do unzip \
$file -d /home/workspace/ccorreia/globin/quality_check/post_alignment/cattle/tmp; \
done

for file in \
`find /home/workspace/ccorreia/globin/quality_check/post_alignment/cattle/tmp \
-name summary.txt`; do more $file >> reports_post-alignment.txt; \
done

for file in \
`find /home/workspace/ccorreia/globin/quality_check/post_alignment/cattle/tmp \
-name fastqc_data.txt`; do head -n 10 $file >> basic_stats_post_alignment.txt; \
done

# Check if all files were processed:
grep -c '##FastQC' basic_stats_post_alignment.txt
grep -c 'Basic Statistics' reports_post-alignment.txt
grep -c 'Analysis complete' fastqc_aligned.sh.00.nohup
grep -c 'Analysis complete' fastqc_aligned.sh.01.nohup

# Transfer all files from tmp to laptop:
scp -r remote_server@rodeo.ucd.ie:/home/workspace/ccorreia/globin/quality_check/post_alignment/cattle/tmp .

# Remove tmp folder from Rodeo:
rm -r tmp

###################################################################
# Summarisation of gene counts with featureCounts for sense genes #
###################################################################

# Required package is featureCounts, which is part of Subread 1.5.1 software,
# consult manual for details:
# http://bioinf.wehi.edu.au/subread-package/SubreadUsersGuide.pdf

# Create working directories:
mkdir -p /home/workspace/ccorreia/globin/Count_summarisation/sense/cattle
cd !$

# Run featureCounts with one sample to check if it is working fine:
featureCounts -a \
/home/workspace/genomes/bostaurus/UMD3.1.1_NCBI/annotation_file/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.gff \
-B -p -C -R -s 1 -T 20 -t gene -g Dbxref -o ./sense_counts.txt \
/home/workspace/ccorreia/globin/STAR-2.5.2b_alignment/cattle/A6511_W-1_F/A6511_W-1_F_Aligned.out.bam

# Create a bash script to run featureCounts on BAM file containing multihits and
# uniquely mapped reads using the stranded parameter:
for file in `find /home/workspace/ccorreia/globin/STAR-2.5.2b_alignment/cattle \
-name *_Aligned.out.bam`; \
do sample=`basename $file | perl -p -e 's/_Aligned.out.bam//'`; \
echo "mkdir /home/workspace/ccorreia/globin/Count_summarisation/sense/cattle/$sample; \
cd /home/workspace/ccorreia/globin/Count_summarisation/sense/cattle/$sample; \
featureCounts -a \
/home/workspace/genomes/bostaurus/UMD3.1.1_NCBI/annotation_file/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.gff \
-B -p -C -R -s 1 -T 20 -t gene -g Dbxref \
-o ${sample}_sense-counts.txt $file" >> sense_count.sh; \
done

# Split and run all scripts on Rodeo:
split -d -l 5 sense_count.sh sense_count.sh.
for script in `ls sense_count.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Check if all files were processed:
grep -c 'Read assignment finished.' sense_count.sh.00.nohup
grep -c 'Read assignment finished.' sense_count.sh.01.nohup

# Create bash script to merge stats info from .featureCounts from all samples
# into a single file:
for file in `find /home/workspace/ccorreia/globin/Count_summarisation/sense/cattle \
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
mkdir /home/workspace/ccorreia/globin/Count_summarisation/sense/cattle/tmp

for file in `find /home/workspace/ccorreia/globin/Count_summarisation/sense/cattle \
-name *sense-counts.txt`; do cp $file \
-t /home/workspace/ccorreia/globin/Count_summarisation/sense/cattle/tmp; \
done

# Transfer all files from tmp to laptop, then remove tmp folder:
scp -r ccorreia@rodeo.ucd.ie:/home/workspace/ccorreia/globin/Count_summarisation/sense/cattle/tmp .

# Remove tmp folder from Rodeo:
rm -r tmp

# Following steps will be performed in R.


































