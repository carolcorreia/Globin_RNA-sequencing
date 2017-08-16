##############################################################################
#      RNA-seq analysis of reticulocyte-derived globin gene transcripts      #
#             (HBA and HBB) in bovine peripheral blood samples               #
#                  --- Linux bioinformatics workflow ---                     #
##############################################################################

# DOI badge: 
# Author: Correia, C.N.
# Version 1.0.0
# Last updated on: 16/08/2017

################################
# Download and files check sum #
################################

# All files were downloaded from MSU in 2013. At the time, they were
# also md5sum checked and renamed.

# File names correspond to:
# AnimalNumber_TimePoint_AnimalGroup_PairedEndTag_LaneNumber_fastq.gz

###########################################
# FastQC quality check of raw FASTQ files #
###########################################

# Required software is FastQC v0.11.5, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and enter the quality check output directory:
mkdir -p $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/pre-filtering
cd !$

# Run FastQC in one file to see if it's working well:
fastqc -o $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/pre-filtering \
--noextract --nogroup -t 2 \
/workspace/storage/kmcloughlin/RNAseqTimeCourse/A6511_W10_P_R1_001.fastq.gz

### Moved this folder to my laptop via SCP
### and checked the HTML report. It worked fine.

# Create a bash script to perform FastQC quality check on all fastq.gz files:
for file in `find /workspace/storage/kmcloughlin/RNAseqTimeCourse/ \
-name *fastq.gz`; do echo "fastqc --noextract --nogroup -t 1 \
-o $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/pre-filtering $file" \
>> fastqc.sh; done;

# Split and run all scripts on Stampede:
split -d -l 70 fastqc.sh fastqc.sh.
for script in `ls fastqc.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Check if all the files were processed:
for file in `ls fastqc.sh.0*.nohup`; \
do more $file | grep "Failed to process file" >> failed_fastqc.txt
done

# Deleted all the HTML files:
rm -r *.html

# Check all output from FastQC:
mkdir $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/pre-filtering/tmp

for file in `ls *_fastqc.zip`; do unzip \
$file -d $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/pre-filtering/tmp; \
done;

for file in \
`find $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/pre-filtering/tmp \
-name summary.txt`; do more $file >> reports_pre-filtering.txt; done

for file in \
`find $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/pre-filtering/tmp \
-name fastqc_data.txt`; do head -n 10 $file >> basic_stats_pre-filtering.txt; \
done

# Remove temporary folder and its files:
rm -rf $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/pre-filtering/tmp


##################################################################
# Adapter-contamination and quality filtering of raw FASTQ files #
##################################################################

# Required software is ngsShoRT (version 2.2). More information can be found
# here: http://research.bioinformatics.udel.edu/genomics/ngsShoRT/index.html

# Create a working directory for filtered reads:
mkdir $HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/
cd $HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/

# Run ngsShoRT in one pair of reads to check if it's working:
nohup perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 20 -mode trim -min_rl 100 \
-pe1 /workspace/storage/kmcloughlin/RNAseqTimeCourse/A6522_W10_U_R1_004.fastq.gz \
-pe2 /workspace/storage/kmcloughlin/RNAseqTimeCourse/A6522_W10_U_R2_004.fastq.gz \
-o /home/ccorreia/scratch/PPDbRNAseqTimeCourse/fastq_sequence/A6522_W10_U_004 \
-methods 5adpt_lqr -5a_f Illumina_PE_adapters.txt -5a_mp 90 -5a_del 0 \
-5a_ins 0 -5a_fmi 100 -5a_axn kr -lqs 20 -lq_p 25 -gzip &

# Create bash scripts to perform filtering of each FASTQ file, keeping the
# sequencing lane information:
for file in `find /workspace/storage/kmcloughlin/RNAseqTimeCourse/ \
-name *R1_001.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/R1(_00.\.fastq.gz)$/R2$1/'`; \
sample=`basename $file | perl -p -e 's/R1(_00.\.fastq.gz)$/001/'`; \
echo "perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 15 -mode trim -min_rl 100 \
-pe1 $file -pe2 $file2 \
-o $HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/$sample \
-methods 5adpt_lqr -5a_f Illumina_PE_adapters.txt -5a_mp 90 \
-5a_del 0 -5a_ins 0 -5a_fmi 100 -5a_axn kr -lqs 20 -lq_p 25 -gzip" \
>> filtering.001.sh; done;

for file in `find /workspace/storage/kmcloughlin/RNAseqTimeCourse/ \
-name *R1_002.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/R1(_00.\.fastq.gz)$/R2$1/'`; \
sample=`basename $file | perl -p -e 's/R1(_00.\.fastq.gz)$/002/'`; \
echo "perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 15 -mode trim -min_rl 100 \
-pe1 $file -pe2 $file2 \
-o $HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/$sample \
-methods 5adpt_lqr -5a_f Illumina_PE_adapters.txt -5a_mp 90 \
-5a_del 0 -5a_ins 0 -5a_fmi 100 -5a_axn kr -lqs 20 -lq_p 25 -gzip" \
>> filtering.002.sh; done;

for file in `find /workspace/storage/kmcloughlin/RNAseqTimeCourse/ \
-name *R1_003.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/R1(_00.\.fastq.gz)$/R2$1/'`; \
sample=`basename $file | perl -p -e 's/R1(_00.\.fastq.gz)$/003/'`; \
echo "perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 15 -mode trim -min_rl 100 \
-pe1 $file -pe2 $file2 \
-o $HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/$sample \
-methods 5adpt_lqr -5a_f Illumina_PE_adapters.txt -5a_mp 90 \
-5a_del 0 -5a_ins 0 -5a_fmi 100 -5a_axn kr -lqs 20 -lq_p 25 -gzip" \
>> filtering.003.sh; done;

for file in `find /workspace/storage/kmcloughlin/RNAseqTimeCourse/ \
-name *R1_004.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/R1(_00.\.fastq.gz)$/R2$1/'`; \
sample=`basename $file | perl -p -e 's/R1(_00.\.fastq.gz)$/004/'`; \
echo "perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 15 -mode trim -min_rl 100 \
-pe1 $file -pe2 $file2 \
-o $HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/$sample \
-methods 5adpt_lqr -5a_f Illumina_PE_adapters.txt -5a_mp 90 \
-5a_del 0 -5a_ins 0 -5a_fmi 100 -5a_axn kr -lqs 20 -lq_p 25 -gzip" \
>> filtering.004.sh; done;

for file in `find /workspace/storage/kmcloughlin/RNAseqTimeCourse/ \
-name *R1_005.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/R1(_00.\.fastq.gz)$/R2$1/'`; \
sample=`basename $file | perl -p -e 's/R1(_00.\.fastq.gz)$/005/'`; \
echo "perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 15 -mode trim -min_rl 100 \
-pe1 $file -pe2 $file2 \
-o $HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/$sample \
-methods 5adpt_lqr -5a_f Illumina_PE_adapters.txt -5a_mp 90 \
-5a_del 0 -5a_ins 0 -5a_fmi 100 -5a_axn kr -lqs 20 -lq_p 25 -gzip" \
>> filtering.005.sh; done;

# Run all scripts on Stampede:
for script in `ls filtering.00*.sh`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Check that all files were processed:
for file in `ls filtering.00*.sh.nohup`; \
do grep -o 'Done-MAIN' $file | wc -l; done

# Compress files with discarded reads:
for file in `find $HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence \
-name extracted_*.txt`; do echo "gzip -9 $file" >> discarded_compression.sh; \
done;

# Split and run all scripts on Stampede:
split -d -l 70 discarded_compression.sh discarded_compression.sh.
for script in `ls discarded_compression.sh.*`
do
chmod 755 $script
nohup ./$script &
done

################################################
# FastQC quality check of filtered FASTQ files #
################################################

# Required software is FastQC v0.11.5, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and enter the quality check output directory:
mkdir $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/post-filtering
cd $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/post-filtering

# Run FastQC in one file to see if it's working well:
fastqc -o $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/post-filtering \
--noextract --nogroup -t 10 \
$HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/A6644_W6_F_001/trimmed_A6644_W6_F_R1_001.fastq.gz

### Moved this folder to my laptop via SCP
### and checked the HTML report. It worked fine.

# Create a bash script to perform FastQC quality check on all filtered
# FASTQ files:
for file in `find $HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/ \
-name *_R*_00*.fastq.gz`; do echo "fastqc --noextract --nogroup -t 1 \
-o $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/post-filtering $file" \
>> fastqc_filt.sh; done;

# Split and run all scripts on Stampede:
split -d -l 75 fastqc_filt.sh fastqc_filt.sh.
for script in `ls fastqc_filt.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Check if all the files were processed:
for file in `ls fastqc_filt.sh.0*.nohup`; \
do more $file | grep "Failed to process file" >> failed_fastqc.txt
done

# Deleted all the HTML files:
rm -r *.html

# Check all output from FastQC:
mkdir $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/post-filtering/tmp

for file in `ls *_fastqc.zip`; do unzip \
$file -d $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/post-filtering/tmp; \
done;

for file in \
`find $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/post-filtering/tmp \
-name summary.txt`; do more $file >> reports_post-filtering.txt; done

for file in \
`find $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/post-filtering/tmp \
-name fastqc_data.txt`; do head -n 10 $file >> basic_stats_post-filtering.txt; \
done

# Remove temporary folder and its files:
rm -rf $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/post-filtering/tmp


### Moved trimmed reads to Rodeo.
### Following steps were conducted in Ubuntu 14.04

#####################################################
# Download reference transcriptome from NCBI RefSeq #
#####################################################

# ASSEMBLY NAME: Bos_taurus_UMD_3.1.1
# ASSEMBLY ACCESSION: GCF_000003055.6

# Create and enter working directory:
mkdir -p /home/workspace/ccorreia/globin/RefSeq/transcriptomes/cattle/source_file
cd !$

# Download and unzip the transcriptome:
nohup wget ftp://ftp.ncbi.nih.gov/genomes/Bos_taurus/RNA/rna.fa.gz &
gunzip rna.fa.gz

# Correct transcript IDs by removing the extra characters 'gi||' and 'ref||':
sed -e 's/gi|.*|ref|//' -e 's/|//' rna.fa > bta_refMrna.fa

##############################################
# Build the transcriptome index using Salmon #
##############################################

# Required software is Salmon 0.8.2, consult manual/tutorial for details:
http://salmon.readthedocs.io/en/latest/

# Enter working directory:
cd /home/workspace/ccorreia/globin/RefSeq/transcriptomes/cattle

# Build an index for quasi-mapping:
nohup salmon index -t \
/home/workspace/ccorreia/globin/RefSeq/transcriptomes/cattle/source_file/bta_refMrna.fa \
-i cattle_index --type quasi -k 31 -p 20 &

#####################################
# Quantify transcripts using Salmon #
#####################################

# Create and enter working directory:
mkdir -p /home/workspace/ccorreia/globin/RefSeq/salmon_quant/cattle
cd !$

# Quantify transcripts from one FASTQ file to check if it works well:
salmon quant -i /home/workspace/ccorreia/globin/RefSeq/transcriptomes/cattle/cattle_index \
-l A --seqBias --gcBias \
-1 \
/home/workspace/ccorreia/globin/fastq_sequence/cattle/trimmed_A6511_W-1_F_R1_001.fastq.gz \
/home/workspace/ccorreia/globin/fastq_sequence/cattle/trimmed_A6511_W-1_F_R1_002.fastq.gz \
/home/workspace/ccorreia/globin/fastq_sequence/cattle/trimmed_A6511_W-1_F_R1_003.fastq.gz \
/home/workspace/ccorreia/globin/fastq_sequence/cattle/trimmed_A6511_W-1_F_R1_004.fastq.gz \
/home/workspace/ccorreia/globin/fastq_sequence/cattle/trimmed_A6511_W-1_F_R1_005.fastq.gz \
-2 \
/home/workspace/ccorreia/globin/fastq_sequence/cattle/trimmed_A6511_W-1_F_R2_001.fastq.gz \
/home/workspace/ccorreia/globin/fastq_sequence/cattle/trimmed_A6511_W-1_F_R2_002.fastq.gz \
/home/workspace/ccorreia/globin/fastq_sequence/cattle/trimmed_A6511_W-1_F_R2_003.fastq.gz \
/home/workspace/ccorreia/globin/fastq_sequence/cattle/trimmed_A6511_W-1_F_R2_004.fastq.gz \
/home/workspace/ccorreia/globin/fastq_sequence/cattle/trimmed_A6511_W-1_F_R2_005.fastq.gz \
-p 15 -o ./A6511_W-1_F 

# Create a bash script to perform quantification of paired FASTQ files sequenced
# over different lanes.
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
echo "salmon quant -i \
/home/workspace/ccorreia/globin/RefSeq/transcriptomes/cattle/cattle_index \
-l A --seqBias --gcBias -1 $file $file2 $file3 $file4 $file5 \
-2 $read1 $read2 $read3 $read4 $read5 \
-p 15 -o ./$sample" \
>> quantify.sh; \
done

# Run script on Rodeo:
chmod 755 quantify.sh
nohup ./quantify.sh > quantify.sh.nohup &

# Append sample name to all quant.sf files:
for file in `find /home/workspace/ccorreia/globin/RefSeq/salmon_quant/cattle \
-name quant.sf`; \
do oldname=`basename $file`; \
newname=`dirname $file | perl -p -e 's/.+(A\d\d\d\d_W\-1_F)/$1/'`; \
path=`dirname $file`; \
mv $file $path/${newname}_$oldname; \
done

# Move all *quant.sf files to a temporary folder:
mkdir /home/workspace/ccorreia/globin/RefSeq/salmon_quant/cattle/cattle_TPM
for file in `find /home/workspace/ccorreia/globin/RefSeq/salmon_quant/cattle \
-name A*quant.sf`; \
do cp $file -t /home/workspace/ccorreia/globin/RefSeq/salmon_quant/cattle/cattle_TPM; \
done

# Transfer all files from Rodeo to laptop:
scp -r ccorreia@servername:/home/workspace/ccorreia/globin/RefSeq/salmon_quant/cattle/cattle_TPM .

# Remove tmp folder from Rodeo:
rm -r cattle_TPM

# Append sample name to all log files:
for file in `find /home/workspace/ccorreia/globin/RefSeq/salmon_quant/cattle/ \
-name salmon_quant.log`; \
do oldname=`basename $file`; \
newname=`dirname $file | perl -p -e 's/.+(A\d\d\d\d_W\-1_F).+/$1/'`; \
path=`dirname $file`; \
mv $file $path/${newname}_$oldname; \
done

# Gather salmon log information from all samples into one file:
for file in `find /home/workspace/ccorreia/globin/RefSeq/salmon_quant/cattle/ \
-name A*salmon_quant.log`; \
do echo echo \
"\`basename $file\` \
\`grep 'likely library type' $file | awk '{print \$12}'\` \
\`grep 'total fragments' $file | awk '{print \$2}'\` \
\`grep 'total reads' $file | awk '{print \$6}'\` \
\`grep 'Mapping rate' $file | awk '{print \$8}'\` >> \
RefSeq_summary_cattle.txt" >> RefSeq_summary_cattle.sh
done

chmod 755 RefSeq_summary_cattle.sh
./RefSeq_summary_cattle.sh

sed -i $'1 i\\\nFile_name Library_type Total_fragments Total_reads Mapping_rate(%)' \
RefSeq_summary_cattle.txt

# Transfer summary file from Rodeo to laptop:
scp -r ccorreia@servername:/home/workspace/ccorreia/globin/RefSeq/salmon_quant/cattle/RefSeq_summary_cattle.txt .

#######################################
# Following steps were performed in R #
#######################################

# Please check this file for next steps: 01-GlobinRNA-seqAnalysis.R









