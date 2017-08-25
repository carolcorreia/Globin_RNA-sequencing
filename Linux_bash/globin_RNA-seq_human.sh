##############################################################################
#      RNA-seq analysis of reticulocyte-derived globin gene transcripts      #
#             (HBA and HBB) in human peripheral blood samples                #
#                  --- Linux bioinformatics workflow ---                     #
##############################################################################

# DOI badge: 
# Author: Correia, C.N.
# Version 1.0.0
# Last updated on: 25/08/2017

#####################################
# Download raw FASTQ files from ENA #
#####################################

# Create and enter working directory:
mkdir $HOME/storage/globin/human_fastq
cd !$

# Get FTP links for FASTQ files from the European Nucleotide Archive (ENA)
# http://www.ebi.ac.uk/ena/data/view/PRJNA232593
# Under the Read Files tab, click on select columns and make sure that the
# only checked option is FASTQ files (FTP). Then, on the download option 
# click on the TEXT link and download the PRJNA232593.txt file.
# Transfer it to Stampede using WinSCP.
# Using notepad++, remove the header on PRJNA232593.txt. Then, search ; and
# replace with \n
# Every FTP link should be in one line.

# Download human data set from ENA (PRJNA232593, Shin et al.2014):
nohup wget -i PRJNA232593.txt &

# Check that all files were downloaded:
grep '.fastq.gz' PRJNA232593.txt | wc -l   # 94 lines
grep '.fastq.gz’ saved' nohup.out | wc -l  # 94 lines

###########################################
# FastQC quality check of raw FASTQ files #
###########################################

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

##################################################################
# Adapter-contamination and quality filtering of raw FASTQ files #
##################################################################

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

# Gather ngsShoRT reports from all samples into one file:
for file in `find $HOME/scratch/globin/fastq_sequence/human \
-name final_PE_report.txt`; \
do echo echo \
"\`dirname $file | perl -p -e 's/.*(SRR\d.+)/\$1/'\` \
\`grep 'Read Pair Count:' $file\` \
\`grep 'Removed PE Pair\* Count:' $file\` >> \
ngsshort_human.txt" >> ngsshort_summary_human.sh
done

chmod 755 ngsshort_summary_human.sh
./ngsshort_summary_human.sh

# Transfer ngsShoRT summary to laptop via SCP:

################################################
# FastQC quality check of filtered FASTQ files #
################################################

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


##################################
# Rename ENA fastq trimmed files #
##################################

### Moved trimmed reads to Rodeo.
### Following steps were conducted in Ubuntu 14.04

# Manually download the "RunInfo Table" file from NCBI, then transfer it
# to Rodeo via WinSCP:
# https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP034732

# Enter working directory in Rodeo:
cd /home/workspace/ccorreia/globin/fastq_sequence/human

# Create a file with two columns: the current folder name and
# the desired new name.
# From the SRA Run file select the following columns:
# 'Run_s', 'treatment_s' internal_id_s, and 'lane_s' 
# and append "Subj" and "L" to tidy the sample information.
awk '{if (NR!=1) {print $5, $5"_"$20"_"$16"_""Subj"$10"_""L"$11}}' \
SraRunTable.txt > sample_info.txt

# Rename the folders with a a custom python script:
python3 /home/workspace/ccorreia/scripts/rename_files.py sample_info.txt

# Create a temporary folder and copy all data to it:
mkdir /home/workspace/ccorreia/globin/fastq_sequence/human/tmp

for file in `find /home/workspace/ccorreia/globin/fastq_sequence/human/SRR* \
-name trimmed_*.fastq.gz`; \
do cp $file \
-t /home/workspace/ccorreia/globin/fastq_sequence/human/tmp; \
done

# Create a second file with two columns: the current FASTQ file name and
# the desired new name.
for file in `find /home/workspace/ccorreia/globin/fastq_sequence/human/SRR* \
-name *_1.fastq.gz`; \
do oldnameR1=`basename $file`; \
oldnameR2=`echo $oldnameR1 | perl -p -e 's/_1\.fastq\.gz/_2\.fastq\.gz/'`; \
pathname=`dirname $file`; \
newname=`basename $pathname | perl -p -e 's/SRR\d+_(.+_.+_.+_L\d)/$1/'`; \
echo "$oldnameR1  ${newname}_1.fastq.gz" >> fastq_names.txt; \
echo "$oldnameR2  ${newname}_2.fastq.gz" >> fastq_names.txt; \
done

# Rename the FASTQ files:
mv fastq_names.txt -t /home/workspace/ccorreia/globin/fastq_sequence/human/tmp
cd /home/workspace/ccorreia/globin/fastq_sequence/human/tmp
python3 /home/workspace/ccorreia/scripts/rename_files.py fastq_names.txt

#####################################################
# Download reference transcriptome from NCBI RefSeq #
#####################################################

# ASSEMBLY NAME: GRCh38.p7
# ASSEMBLY ACCESSION: GCF_000001405.33

# Create and enter working directory:
mkdir -p /home/workspace/ccorreia/globin/RefSeq/transcriptomes/human/source_file
cd !$

# Download and unzip the transcriptome:
nohup wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/RNA/rna.fa.gz &
gunzip rna.fa.gz

# Correct transcript IDs by removing the extra characters 'gi||' and 'ref||':
sed -e 's/gi|.*|ref|//' -e 's/|//' rna.fa > hsa_refMrna.fa

##############################################
# Build the transcriptome index using Salmon #
##############################################

# Required software is Salmon 0.8.2, consult manual/tutorial for details:
http://salmon.readthedocs.io/en/latest/

# Enter working directory:
cd /home/workspace/ccorreia/globin/RefSeq/transcriptomes/human

# Build an index for quasi-mapping:
nohup salmon index -t \
/home/workspace/ccorreia/globin/RefSeq/transcriptomes/human/source_file/hsa_refMrna.fa \
-i human_index --type quasi -k 31 -p 20 &

#####################################
# Quantify transcripts using Salmon #
#####################################

# Create and enter working directory:
mkdir -p /home/workspace/ccorreia/globin/RefSeq/salmon_quant/human
cd !$

# Quantify transcripts from one FASTQ file to check if it works well:
salmon quant -i /home/workspace/ccorreia/globin/RefSeq/transcriptomes/human/human_index \
-l A --seqBias --gcBias \
-1 /home/workspace/ccorreia/globin/fastq_sequence/human/renamed_fastq/GD_TRUE_Subj12_L1_1.fastq.gz \
-2 /home/workspace/ccorreia/globin/fastq_sequence/human/renamed_fastq/GD_TRUE_Subj12_L1_2.fastq.gz \
-p 15 -o ./GD_FALSE_Subj12

# Create a bash script to perform quantification of paired FASTQ files sequenced
# over different lanes.
# Lanes 1 and 4:
for file in `find /home/workspace/ccorreia/globin/fastq_sequence/human/renamed_fastq \
-name *L1_1.fastq.gz`; \
do read1=`echo $file | perl -p -e 's/_L1_1\.fastq\.gz/\_L1_2\.fastq\.gz/'`; \
file2=`echo $file | perl -p -e 's/\_L1_1\.fastq\.gz/\_L4_1\.fastq\.gz/'`; \
read2=`echo $file2 | perl -p -e 's/\_L4_1\.fastq\.gz/\_L4_2\.fastq\.gz/'`; \
sample=`basename $file | perl -p -e 's/(.+_.+_Subj\d+)_L\d_\d\.fastq\.gz/$1/'`; \
echo "salmon quant -i /home/workspace/ccorreia/globin/RefSeq/transcriptomes/human/human_index \
-l A --seqBias --gcBias -1 $file $file2 -2 $read1 $read2 \
-p 15 -o ./$sample" \
>> quant.sh; \
done

# Delete line for Subj 12. This was the only sample not
# sequenced over two lanes:
sed '/Subj12/d' ./quant.sh > quantify.sh
rm quant.sh

# Lanes 2 and 5:
for file in `find /home/workspace/ccorreia/globin/fastq_sequence/human/renamed_fastq \
-name *L2_1.fastq.gz`; \
do read1=`echo $file | perl -p -e 's/_L2_1\.fastq\.gz/\_L2_2\.fastq\.gz/'`; \
file2=`echo $file | perl -p -e 's/\_L2_1\.fastq\.gz/\_L5_1\.fastq\.gz/'`; \
read2=`echo $file2 | perl -p -e 's/\_L5_1\.fastq\.gz/\_L5_2\.fastq\.gz/'`; \
sample=`basename $file | perl -p -e 's/(.+_.+_Subj\d+)_L\d_\d\.fastq\.gz/$1/'`; \
echo "salmon quant -i /home/workspace/ccorreia/globin/RefSeq/transcriptomes/human/human_index \
-l A --seqBias --gcBias -1 $file $file2 -2 $read1 $read2 \
-p 15 -o ./$sample" \
>> quantify.sh; \
done

# Split and run all scripts on Rodeo:
split -d -l 12 quantify.sh quantify.sh.
for script in `ls quantify.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Append sample name to all quant.sf files to temporary folder:
for file in `find /home/workspace/ccorreia/globin/RefSeq/salmon_quant/human \
-name quant.sf`; \
do oldname=`basename $file`; \
newname=`dirname $file | perl -p -e 's/.+(GD|NGD_.+_Subj\d+)/$1/'`; \
path=`dirname $file`; \
mv $file $path/${newname}_$oldname; \
done

# Move all *quant.sf files to a temporary folder:
mkdir /home/workspace/ccorreia/globin/RefSeq/salmon_quant/human/human_TPM
for file in `find /home/workspace/ccorreia/globin/RefSeq/salmon_quant/human \
-name *_Subj*quant.sf`; \
do cp $file -t /home/workspace/ccorreia/globin/RefSeq/salmon_quant/human/human_TPM; \
done

# Transfer all files from Rodeo to laptop:
scp -r ccorreia@remoteserver:/home/workspace/ccorreia/globin/RefSeq/salmon_quant/human/human_TPM .

# Remove tmp folder from Rodeo:
rm -r human_TPM

# Append sample name to all log files:
for file in `find /home/workspace/ccorreia/globin/RefSeq/salmon_quant/human/ \
-name salmon_quant.log`; \
do oldname=`basename $file`; \
newname=`dirname $file | perl -p -e 's/.+(.+D_.+_Subj\d+).+/$1/'`; \
path=`dirname $file`; \
mv $file $path/${newname}_$oldname; \
done

# Gather salmon log information from all samples into one file:
for file in `find /home/workspace/ccorreia/globin/RefSeq/salmon_quant/human/ \
-name *D_*salmon_quant.log`; \
do echo echo \
"\`basename $file\` \
\`grep 'likely library type' $file | awk '{print \$12}'\` \
\`grep 'total fragments' $file | awk '{print \$2}'\` \
\`grep 'total reads' $file | awk '{print \$6}'\` \
\`grep 'Mapping rate' $file | awk '{print \$8}'\` >> \
RefSeq_summary_human.txt" >> RefSeq_summary_human.sh
done

chmod 755 RefSeq_summary_human.sh
./RefSeq_summary_human.sh

sed -i $'1 i\\\nFile_name Library_type Total_fragments Total_reads Mapping_rate(%)' \
RefSeq_summary_human.txt

# Transfer summary file from Rodeo to laptop:
scp ccorreia@remoteserver:/home/workspace/ccorreia/globin/RefSeq/salmon_quant/human/RefSeq_summary_human.txt .

#######################################
# Following steps were performed in R #
#######################################

# Please check this file for next steps: 01-GlobinRNA-seqAnalysis.R


