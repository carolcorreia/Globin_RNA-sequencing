##############################################################################
#      RNA-seq analysis of reticulocyte-derived globin gene transcripts      #
#             (HBA and HBB) in equine peripheral blood samples               #
#                           equCab3 - NCBI - 2018                            #
#                  --- Linux bioinformatics workflow ---                     #
##############################################################################

# DOI badge: 
# Author: Correia, C.N.
# Version 1.0.0
# Last updated on: 25/05/2018

#####################################
# Download raw FASTQ files from ENA #
#####################################

# See 4-globin_RNA-seq_horse_UCSC.sh

###########################################
# FastQC quality check of raw FASTQ files #
###########################################

# See 4-globin_RNA-seq_horse_UCSC.sh

##################################################################
# Adapter-contamination and quality filtering of raw FASTQ files #
##################################################################

# See 4-globin_RNA-seq_horse_UCSC.sh

################################################
# FastQC quality check of filtered FASTQ files #
################################################

# See 4-globin_RNA-seq_horse_UCSC.sh

##############################################
# Download reference transcriptome from NCBI #
##############################################

# ASSEMBLY NAME: EquCab3.0
# ASSEMBLY ACCESSION: GCF_002863925.1

# From NCBI Equus caballus Annotation Release 103
# https://www.ncbi.nlm.nih.gov/genome/annotation_euk/Equus_caballus/103/

# Create and enter working directory:
mkdir -p /home/workspace/ccorreia/globin/RefSeq/transcriptomes/horse_2018/source_file
cd !$

# Download and unzip the transcriptome:
nohup wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Equus_caballus/RNA/rna.fa.gz &
gunzip rna.fa.gz

# Correct transcript IDs by removing the extra characters 'ref||':
sed -e 's/ref|//' -e 's/|//' rna.fa > eca_refMrna2018.fa

##############################################
# Build the transcriptome index using Salmon #
##############################################

# Required software is Salmon 0.8.2, consult manual/tutorial for details:
http://salmon.readthedocs.io/en/latest/

# Enter working directory:
cd /home/workspace/ccorreia/globin/RefSeq/transcriptomes/horse_2018

# Build an index for quasi-mapping:
nohup salmon index -t \
/home/workspace/ccorreia/globin/RefSeq/transcriptomes/horse_2018/source_file/eca_refMrna2018.fa \
-i horse_index --type quasi -k 31 -p 20 &


#####################################
# Quantify transcripts using Salmon #
#####################################

# Create and enter working directory:
mkdir -p /home/workspace/ccorreia/globin/RefSeq/salmon_quant/horse_2018
cd !$

# Quantify transcripts from one FASTQ file to check if it works well:
salmon quant -i /home/workspace/ccorreia/globin/RefSeq/transcriptomes/horse_2018/horse_index \
-l A  --seqBias --gcBias \
-r /home/workspace/ccorreia/globin/fastq_sequence/horse/SRR3671009/trimmed_SRR3671009.fastq.gz \
-p 15 -o ./SRR3671009

# Create a bash script to perform quantification of FASTQ files sequenced:
for file in `find /home/workspace/ccorreia/globin/fastq_sequence/horse \
-name *.fastq.gz`; \
do sample=`basename $file | perl -p -e 's/.fastq.gz//'`; \
echo "salmon quant -i /home/workspace/ccorreia/globin/RefSeq/transcriptomes/horse_2018/horse_index \
-l A --seqBias --gcBias -r $file -p 15 -o ./$sample" \
>> quantify.sh; \
done

# Run script on Rodeo:
chmod 755 quantify.sh
nohup ./quantify.sh > quantify.sh.nohup &

# Append sample name to all quant.sf files to temporary folder:
for file in `find /home/workspace/ccorreia/globin/RefSeq/salmon_quant/horse_2018 \
-name quant.sf`; \
do sample=`dirname $file | perl -p -e 's/.+_(SRR\d+)/$1/'`; \
oldname=`basename $file`; \
path=`dirname $file`; \
mv $file $path/${sample}_$oldname; \
done

# Move all *quant.sf files to a temporary folder:
mkdir /home/workspace/ccorreia/globin/RefSeq/salmon_quant/horse_2018/horse_TPM
for file in `find /home/workspace/ccorreia/globin/RefSeq/salmon_quant/horse_2018 \
-name SRR*quant.sf`; \
do cp $file -t /home/workspace/ccorreia/globin/RefSeq/salmon_quant/horse_2018/horse_TPM; \
done

# Transfer all files from Rodeo to laptop:
scp -r ccorreia@remoteserver:/home/workspace/ccorreia/globin/RefSeq/salmon_quant/horse_2018/horse_TPM .

# Remove tmp folder from Rodeo:
rm -r horse_TPM

# Append sample name to all log files:
for file in `find /home/workspace/ccorreia/globin/RefSeq/salmon_quant/horse_2018 \
-name salmon_quant.log`; \
do sample=`dirname $file | perl -p -e 's/.+_(SRR\d+).+/$1/'`; \
oldname=`basename $file`; \
path=`dirname $file`; \
mv $file $path/${sample}_$oldname; \
done

# Gather salmon log information from all samples into one file:
for file in `find /home/workspace/ccorreia/globin/RefSeq/salmon_quant/horse_2018 \
-name SRR*salmon_quant.log`; \
do echo echo \
"\`basename $file\` \
\`grep 'likely library type' $file | awk '{print \$12}'\` \
\`grep 'total fragments' $file | awk '{print \$2}'\` \
\`grep 'total reads' $file | awk '{print \$6}'\` \
\`grep 'Mapping rate' $file | awk '{print \$8}'\` >> \
NCBI2018_summary_horse.txt" >> NCBI2018_summary_horse.sh
done

chmod 755 NCBI2018_summary_horse.sh
./NCBI2018_summary_horse.sh

sed -i $'1 i\\\nFile_name Library_type Total_fragments Total_reads Mapping_rate(%)' \
NCBI2018_summary_horse.txt

# Transfer summary file from Rodeo to laptop:
scp ccorreia@remoteserver:/home/workspace/ccorreia/globin/RefSeq/salmon_quant/horse_2018/NCBI2018_summary_horse.txt .

#######################################
# Following steps were performed in R #
#######################################

# Please check this file sfor following steps: 09-Globin-2018_analysis-Part_A.R







