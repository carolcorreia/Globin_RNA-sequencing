# Extra analysis for Equine data with NCBI RefSeq

# DOI badge: 
# Author: Correia, C.N.
# Version 1.0.0
# Last updated on: 23/08/2017

#####################################################
# Download reference transcriptome from NCBI RefSeq #
#####################################################

# Sep. 2007 assembly of the horse genome (equCab2, Broad Institute EquCab2).
# http://hgdownload.soe.ucsc.edu/goldenPath/equCab2/bigZips/

# Create and enter working directory:
mkdir -p /home/workspace/ccorreia/globin/horse_extra/transcriptome/source_file
cd !$

# Download and unzip the transcriptome:
nohup wget ftp://ftp.ncbi.nih.gov/genomes/Equus_caballus/RNA/rna.fa.gz &
gunzip rna.fa.gz

# Correct transcript IDs by removing the extra characters 'gi||' and 'ref||':
sed -e 's/gi|.*|ref|//' -e 's/|//' rna.fa > eca_refMrna.fa

##############################################
# Build the transcriptome index using Salmon #
##############################################

# ASSEMBLY NAME: EquCab2.0
# ASSEMBLY ACCESSION: GCF_000002305.2

# Required software is Salmon 0.8.2, consult manual/tutorial for details:
http://salmon.readthedocs.io/en/latest/

# Enter working directory:
cd /home/workspace/ccorreia/globin/horse_extra/transcriptome

# Build an index for quasi-mapping:
nohup salmon index -t \
/home/workspace/ccorreia/globin/horse_extra/transcriptome/source_file/eca_refMrna.fa \
-i horse_index --type quasi -k 31 -p 20 &

#####################################
# Quantify transcripts using Salmon #
#####################################

# Create and enter working directory:
mkdir -p /home/workspace/ccorreia/globin/horse_extra/salmon_quant
cd !$

# Quantify transcripts from one FASTQ file to check if it works well:
salmon quant -i /home/workspace/ccorreia/globin/horse_extra/transcriptome/horse_index \
-l A  --seqBias --gcBias \
-r /home/workspace/ccorreia/globin/fastq_sequence/horse/SRR3671009/trimmed_SRR3671009.fastq.gz \
-p 15 -o ./SRR3671009

# Create a bash script to perform quantification of FASTQ files sequenced:
for file in `find /home/workspace/ccorreia/globin/fastq_sequence/horse \
-name *.fastq.gz`; \
do sample=`basename $file | perl -p -e 's/.fastq.gz//'`; \
echo "salmon quant -i /home/workspace/ccorreia/globin/horse_extra/transcriptome/horse_index \
-l A --seqBias --gcBias -r $file -p 20 -o ./$sample" \
>> quantify.sh; \
done

# Run script on Rodeo:
split -d -l 19 quantify.sh quantify.sh.
for script in `ls quantify.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Append sample name to all quant.sf files to temporary folder:
for file in `find /home/workspace/ccorreia/globin/horse_extra/salmon_quant \
-name quant.sf`; \
do sample=`dirname $file | perl -p -e 's/.+_(SRR\d+)/$1/'`; \
oldname=`basename $file`; \
path=`dirname $file`; \
mv $file $path/${sample}_$oldname; \
done

# Move all *quant.sf files to a temporary folder:
mkdir /home/workspace/ccorreia/globin/horse_extra/salmon_quant/horse_TPM
for file in `find /home/workspace/ccorreia/globin/horse_extra/salmon_quant \
-name SRR*quant.sf`; \
do cp $file -t /home/workspace/ccorreia/globin/horse_extra/salmon_quant/horse_TPM; \
done

# Transfer all files from Rodeo to laptop:
scp -r ccorreia@remoteserver:/home/workspace/ccorreia/globin/horse_extra/salmon_quant/horse_TPM .

# Remove tmp folder from Rodeo:
rm -r horse_TPM

# Append sample name to all log files:
for file in `find /home/workspace/ccorreia/globin/horse_extra/salmon_quant \
-name salmon_quant.log`; \
do sample=`dirname $file | perl -p -e 's/.+_(SRR\d+).+/$1/'`; \
oldname=`basename $file`; \
path=`dirname $file`; \
mv $file $path/${sample}_$oldname; \
done

# Gather salmon log information from all samples into one file:
for file in `find /home/workspace/ccorreia/globin/horse_extra/salmon_quant \
-name SRR*salmon_quant.log`; \
do echo echo \
"\`basename $file\` \
\`grep 'likely library type' $file | awk '{print \$12}'\` \
\`grep 'total fragments' $file | awk '{print \$2}'\` \
\`grep 'total reads' $file | awk '{print \$6}'\` \
\`grep 'Mapping rate' $file | awk '{print \$8}'\` >> \
Extra_summary_horse.txt" >> Extra_summary_horse.sh
done

chmod 755 Extra_summary_horse.sh
./Extra_summary_horse.sh

sed -i $'1 i\\\nFile_name Library_type Total_fragments Total_reads Mapping_rate(%)' \
Extra_summary_horse.txt

# Transfer summary file from Rodeo to laptop:
scp ccorreia@remoteserver:/home/workspace/ccorreia/globin/horse_extra/salmon_quant/Extra_summary_horse.txt .

#######################################
# Following steps were performed in R #
#######################################

# Please check this file for following steps: horse_extra_analysis.R







