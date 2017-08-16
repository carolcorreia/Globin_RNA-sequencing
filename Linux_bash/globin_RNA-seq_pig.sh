##############################################################################
#      RNA-seq analysis of reticulocyte-derived globin gene transcripts      #
#            (HBA and HBB) in porcine peripheral blood samples               #
#                  --- Linux bioinformatics workflow ---                     #
##############################################################################

# DOI badge: 
# Author: Correia, C.N.
# Version 1.0.0
# Last updated on: 16/08/2017

############################
# Download raw FASTQ files #
############################

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


###########################################
# FastQC quality check of raw FASTQ files #
###########################################

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


###############################
# Trimming of raw FASTQ files #
###############################

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


###############################################
# FastQC quality check of trimmed FASTQ files #
###############################################

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


### Moved trimmed reads to Rodeo.
### Following steps were conducted in Ubuntu 14.04


#####################################################
# Download reference transcriptome from NCBI RefSeq #
#####################################################

# ASSEMBLY NAME: Sscrofa11.1
# ASSEMBLY ACCESSION: GCF_000003025.6

# Create and enter working directory:
mkdir -p /home/workspace/ccorreia/globin/RefSeq/transcriptomes/pig/source_file
cd !$

# Download and unzip the transcriptome:
nohup wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Sus_scrofa/RNA/rna.fa.gz &
gunzip rna.fa.gz

# Correct transcript IDs by removing the extra characters 'ref||':
sed -e 's/ref//' -e 's/|//g' rna.fa > ssc_refrna.fa

##############################################
# Build the transcriptome index using Salmon #
##############################################

# Required software is Salmon 0.8.2, consult manual/tutorial for details:
http://salmon.readthedocs.io/en/latest/

# Enter working directory:
cd /home/workspace/ccorreia/globin/RefSeq/transcriptomes/pig

# Build an index for quasi-mapping:
nohup salmon index -t \
/home/workspace/ccorreia/globin/RefSeq/transcriptomes/pig/source_file/ssc_refrna.fa \
-i pig_index --type quasi -k 31 -p 20 &

#####################################
# Quantify transcripts using Salmon #
#####################################

# Create and enter working directory:
mkdir -p /home/workspace/ccorreia/globin/RefSeq/salmon_quant/pig
cd !$

# Quantify transcripts from one FASTQ file to check if it works well:
salmon quant -i /home/workspace/ccorreia/globin/RefSeq/transcriptomes/pig/pig_index \
-l A --seqBias --gcBias \
-1 /home/workspace/ccorreia/globin/fastq_sequence/pig/HI.0751.004.Index_12.GCswine-5037-28DPI-WB-7413C/trimmed_HI.0751.004.Index_12.GCswine-5037-28DPI-WB-7413C-mRNA_R1.fastq.gz \
/home/workspace/ccorreia/globin/fastq_sequence/pig/HI.0751.005.Index_12.GCswine-5037-28DPI-WB-7413C/trimmed_HI.0751.005.Index_12.GCswine-5037-28DPI-WB-7413C-mRNA_R1.fastq.gz \
-2 /home/workspace/ccorreia/globin/fastq_sequence/pig/HI.0751.004.Index_12.GCswine-5037-28DPI-WB-7413C/trimmed_HI.0751.004.Index_12.GCswine-5037-28DPI-WB-7413C-mRNA_R2.fastq.gz \
/home/workspace/ccorreia/globin/fastq_sequence/pig/HI.0751.005.Index_12.GCswine-5037-28DPI-WB-7413C/trimmed_HI.0751.005.Index_12.GCswine-5037-28DPI-WB-7413C-mRNA_R2.fastq.gz \
-p 15 -o ./7413C

# Create a bash script to perform quantification of paired FASTQ files sequenced
# over different lanes.
# Lanes 002 or 003 (contain only one pair for each library):
for file in `find /home/workspace/ccorreia/globin/fastq_sequence/pig \
-name trimmed_*.00[23].*_R1.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/\_R1\.fastq\.gz/\_R2\.fastq\.gz/'`; \
sample=`basename $file | perl -p -e 's/trimmed_.*(7\d\d\d.+)\-mRNA_R1\.fastq\.gz/$1/'`; \
echo "salmon quant -i /home/workspace/ccorreia/globin/RefSeq/transcriptomes/pig/pig_index \
-l A --seqBias --gcBias -1 $file -2 $file2 \
-p 15 -o ./$sample" \
>> quantify.sh; \
done

# Lanes 004 and 005 (contain duplicate pairs for each library):
for file in `find /home/workspace/ccorreia/globin/fastq_sequence/pig \
-name trimmed_*.004.*_R1.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/(\.004\.)/\.005\./g'`; \
read1=`echo $file | perl -p -e 's/\_R1\.fastq\.gz/\_R2\.fastq\.gz/'`; \
read2=`echo $file2 | perl -p -e 's/\_R1\.fastq\.gz/\_R2\.fastq\.gz/'`; \
sample=`basename $file | perl -p -e 's/trimmed_.*(7\d\d\d.+)\-mRNA_R1\.fastq\.gz/$1/'`; \
echo "salmon quant -i /home/workspace/ccorreia/globin/RefSeq/transcriptomes/pig/pig_index \
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
for file in `find /home/workspace/ccorreia/globin/RefSeq/salmon_quant/pig \
-name quant.sf`; \
do oldname=`basename $file`; \
newname=`dirname $file | perl -p -e 's/.+(\d\d\d\d.+)/$1/'`; \
path=`dirname $file`; \
mv $file $path/${newname}_$oldname; \
done

# Move all *quant.sf files to a temporary folder:
mkdir /home/workspace/ccorreia/globin/RefSeq/salmon_quant/pig/pig_TPM
for file in `find /home/workspace/ccorreia/globin/RefSeq/salmon_quant/pig \
-name 7*quant.sf`; \
do cp $file -t /home/workspace/ccorreia/globin/RefSeq/salmon_quant/pig/pig_TPM; \
done

# Transfer all files from Rodeo to laptop:
scp -r ccorreia@remoteserver:/home/workspace/ccorreia/globin/RefSeq/salmon_quant/pig/pig_TPM .

# Remove tmp folder from Rodeo:
rm -r pig_TPM

# Append sample name to all log files:
for file in `find /home/workspace/ccorreia/globin/RefSeq/salmon_quant/pig/ \
-name salmon_quant.log`; \
do oldname=`basename $file`; \
newname=`dirname $file | perl -p -e 's/.+(7\d\d\d\w+).+/$1/'`; \
path=`dirname $file`; \
mv $file $path/${newname}_$oldname; \
done

# Gather salmon log information from all samples into one file:
for file in `find /home/workspace/ccorreia/globin/RefSeq/salmon_quant/pig/ \
-name 7*salmon_quant.log`; \
do echo echo \
"\`basename $file\` \
\`grep 'likely library type' $file | awk '{print \$12}'\` \
\`grep 'total fragments' $file | awk '{print \$2}'\` \
\`grep 'total reads' $file | awk '{print \$6}'\` \
\`grep 'Mapping rate' $file | awk '{print \$8}'\` >> \
RefSeq_summary_pig.txt" >> RefSeq_summary_pig.sh
done

chmod 755 RefSeq_summary_pig.sh
./RefSeq_summary_pig.sh

sed -i $'1 i\\\nFile_name Library_type Total_fragments Total_reads Mapping_rate(%)' \
RefSeq_summary_pig.txt

# Transfer summary file from Rodeo to laptop:
scp -r ccorreia@remoteserver:/home/workspace/ccorreia/globin/RefSeq/salmon_quant/pig/RefSeq_summary_pig.txt .

#######################################
# Following steps were performed in R #
#######################################

# Please check this file for next steps: 01-GlobinRNA-seqAnalysis.R


