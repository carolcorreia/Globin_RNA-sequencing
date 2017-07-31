##############################################################################
#      RNA-seq analysis of reticulocyte-derived globin gene transcripts      #
#             (HBA and HBB) in horse peripheral blood samples                #
#                  --- Linux bioinformatics workflow ---                     #
##############################################################################

# DOI badge: 
# Author: Correia, C.N.
# Version 1.0.0
# Last updated on: 25/07/2017

#####################################
# Download raw FASTQ files from ENA #
#####################################

# Create and enter working directory:
mkdir -p /home/workspace/ccorreia/globin/raw/horse
cd !$

# Get ftp links for FASTQ files from the European Nucleotide Archive (ENA)
# http://www.ebi.ac.uk/ena/data/view/PRJNA325820
# Under the Read Files tab, click on select columns and make sure that the
# only checked option is FASTQ files (FTP). Then, on the download option 
# click on the TEXT link and download the PRJNA325820.txt file.
# Using Sublime Text, remove the header on PRJNA325820.txt file.
# Every FTP link should be in one line.
# Transfer it to Rodeo via SCP.

# Download horse data set from ENA (PRJNA325820, Ropka-Molik et al.2017):
nohup wget -i PRJNA325820.txt &

# Check that all files were downloaded:
grep '.fastq.gz' PRJNA325820.txt | wc -l   # 37 lines
grep '.fastq.gz’ saved' nohup.out | wc -l  # 37 lines

###########################################
# FastQC quality check of raw FASTQ files #
###########################################

# Required software is FastQC v0.11.5, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and enter the quality check output directory:
mkdir -p /home/workspace/ccorreia/globin/quality_check/pre-filtering/horse
cd !$

# Run FastQC in one file to see if it's working well:
fastqc -o /home/workspace/ccorreia/globin/quality_check/pre-filtering/horse \
--noextract --nogroup -t 5 \
/home/workspace/ccorreia/globin/raw/horse/SRR3671009.fastq.gz

### Moved this folder to my laptop via SCP
### and checked the HTML report. It worked fine.

# Create a bash script to perform FastQC quality check on all fastq.gz files:
for file in `find /home/workspace/ccorreia/globin/raw/horse \
-name *fastq.gz`; do echo "fastqc --noextract --nogroup -t 10 \
-o /home/workspace/ccorreia/globin/quality_check/pre-filtering/horse $file" \
>> fastqc.sh; 
done

# Run script on Rodeo:
chmod 755 fastqc.sh
nohup ./fastqc.sh > fastqc.sh.nohup &

# Deleted all the HTML files:
rm -r *.html

# Check if all the files were processed:
ls -l | grep fastqc.zip | wc -l
more fastqc.sh.nohup | grep "Failed to process file" >> failed_fastqc.txt

### Copied all .zip files to my laptop via SCP
### and checked the HTML reports. It worked fine.

# Check all output from FastQC:
mkdir /home/workspace/ccorreia/globin/quality_check/pre-filtering/horse/tmp

for file in `ls *_fastqc.zip`; do unzip \
$file -d /home/workspace/ccorreia/globin/quality_check/pre-filtering/horse/tmp; \
done

for file in \
`find /home/workspace/ccorreia/globin/quality_check/pre-filtering/horse/tmp \
-name summary.txt`; do more $file >> reports_pre-filtering.txt; \
done

grep 'Adapter Content' reports_pre-filtering.txt >> adapter_content.txt
wc -l adapter_content.txt

grep FAIL adapter_content.txt | wc -l

for file in \
`find /home/workspace/ccorreia/globin/quality_check/pre-filtering/horse/tmp \
-name fastqc_data.txt`; do head -n 10 $file >> basic_stats_pre-filtering.txt; \
done

# Remove temporary folder and its files:
rm -r /home/workspace/ccorreia/globin/quality_check/pre-filtering/horse/tmp

##################################################################
# Adapter-contamination and quality filtering of raw FASTQ files #
##################################################################

# Required software is ngsShoRT (version 2.2). More information can be found
# here: http://research.bioinformatics.udel.edu/genomics/ngsShoRT/index.html

# Create a working directory for filtered reads:
mkdir /home/workspace/ccorreia/globin/fastq_sequence/horse
cd !$

# Copy Illumina adpaters file into working directory:
cp /home/workspace/ccorreia/globin/Illumina_PE_adapters.txt .

# Run ngsShoRT in one FASTQ file to check if it's working:
nohup perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 5 -mode trim -min_rl 70 \
-se /home/workspace/ccorreia/globin/raw/horse/SRR3671009.fastq.gz \
-o /home/workspace/ccorreia/globin/fastq_sequence/horse \
-methods 5adpt_lqr -5a_f Illumina_PE_adapters.txt -5a_mp 90 -5a_del 0 \
-5a_ins 0 -5a_fmi 100 -5a_axn kr -lqs 20 -lq_p 25 -gzip &

# Create bash script to perform filtering of each FASTQ file:
for file in `find /home/workspace/ccorreia/globin/raw/horse \
-name *.fastq.gz`; \
do sample=`basename $file | perl -p -e 's/.fastq.gz//'`; \
echo "perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 5 -mode trim -min_rl 70 \
-se $file \
-o /home/workspace/ccorreia/globin/fastq_sequence/horse/$sample \
-methods 5adpt_lqr -5a_f Illumina_PE_adapters.txt -5a_mp 90 -5a_del 0 \
-5a_ins 0 -5a_fmi 100 -5a_axn kr -lqs 20 -lq_p 25 -gzip" \
>> filtering.sh; \
done

# Split and run all scripts on Rodeo:
split -d -l 20 filtering.sh filtering.sh.
for script in `ls filtering.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Check that all folders were created:
ls -l | grep SRR | wc -l

# Check that all files were processed:
for file in `ls filtering.sh.0*.nohup`; \
do grep -o 'Done-MAIN' $file | wc -l; done

# Compress files of removed reads:
for file in `find /home/workspace/ccorreia/globin/fastq_sequence/horse \
-name extracted*.txt`; \
do echo "gzip -9 $file" >> discarded_compression.sh; \
done

# Run script on Rodeo:
chmod 755 discarded_compression.sh
nohup ./discarded_compression.sh &


################################################
# FastQC quality check of filtered FASTQ files #
################################################

# Required software is FastQC v0.11.5, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and enter the quality check output directory:
mkdir -p /home/workspace/ccorreia/globin/quality_check/post-filtering/horse
cd !$

# Run FastQC in one file to see if it's working well:
fastqc -o /home/workspace/ccorreia/globin/quality_check/post-filtering/horse \
--noextract --nogroup -t 5 \
/home/workspace/ccorreia/globin/fastq_sequence/horse/SRR3671009/trimmed_SRR3671009.fastq.gz

### Moved this folder to my laptop via SCP
### and checked the HTML report. It worked fine.

# Create a bash script to perform FastQC quality check on all fastq.gz files:
for file in `find /home/workspace/ccorreia/globin/fastq_sequence/horse \
-name trimmed_*fastq.gz`; do echo "fastqc --noextract --nogroup -t 10 \
-o /home/workspace/ccorreia/globin/quality_check/post-filtering/horse $file" \
>> fastqc.sh; 
done

# Run script on Rodeo:
chmod 755 fastqc.sh
nohup ./fastqc.sh > fastqc.sh.nohup &

# Deleted all the HTML files:
rm -r *.html

# Check if all the files were processed:
ls -l | grep fastqc.zip | wc -l
more fastqc.sh.nohup | grep "Failed to process file" >> failed_fastqc.txt

### Copied all .zip files to my laptop via SCP
### and checked the HTML reports. It worked fine.

# Check all output from FastQC:
mkdir /home/workspace/ccorreia/globin/quality_check/pre-filtering/horse/tmp

for file in `ls *_fastqc.zip`; do unzip \
$file -d /home/workspace/ccorreia/globin/quality_check/pre-filtering/horse/tmp; \
done

for file in \
`find /home/workspace/ccorreia/globin/quality_check/pre-filtering/horse/tmp \
-name summary.txt`; do more $file >> reports_pre-filtering.txt; \
done

grep 'Adapter Content' reports_pre-filtering.txt >> adapter_content.txt
wc -l adapter_content.txt

grep FAIL adapter_content.txt | wc -l

for file in \
`find /home/workspace/ccorreia/globin/quality_check/pre-filtering/horse/tmp \
-name fastqc_data.txt`; do head -n 10 $file >> basic_stats_pre-filtering.txt; \
done

# Remove temporary folder and its files:
rm -r /home/workspace/ccorreia/globin/quality_check/pre-filtering/horse/tmp






















############################################################
# Download reference transcriptome from Ensembl release 88 #
############################################################

# Create and enter working directory:
mkdir -p /home/workspace/ccorreia/globin/transcriptomes/horse/source_file
cd !$

# Download and unzip the transcriptome:
nohup wget ftp://ftp.ensembl.org/pub/release-88/fasta/equus_caballus/cdna/Equus_caballus.EquCab2.cdna.all.fa.gz &
gunzip Equus_caballus.EquCab2.cdna.all.fa.gz

