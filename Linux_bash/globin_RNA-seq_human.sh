##############################################################################
#      RNA-seq analysis of reticulocyte-derived globin gene transcripts      #
#             (HBA and HBB) in human peripheral blood samples                #
#                  --- Linux bioinformatics workflow ---                     #
##############################################################################

# DOI badge: 
# Author: Correia, C.N.
# Version 1.0.0
# Last updated on: 19/05/2017

############################################################
# Download reference transcriptome from Ensembl release 88 #
############################################################

# Create and enter working directory:
mkdir -p /home/workspace/ccorreia/globin/transcriptomes/human/source_file
cd !$

# Download and unzip the transcriptome:
nohup wget ftp://ftp.ensembl.org/pub/release-88/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz &
gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz


##############################################
# Build the transcriptome index using Salmon #
##############################################

# Required software is Salmon 0.8.2, consult manual/tutorial for details:
http://salmon.readthedocs.io/en/latest/

# Enter working directory:
cd /home/workspace/ccorreia/globin/transcriptomes/human

# Build an index for quasi-mapping:
nohup salmon index -t \
/home/workspace/ccorreia/globin/transcriptomes/human/source_file/Homo_sapiens.GRCh38.cdna.all.fa \
-i human_index --type quasi -k 31 -p 20 &


#####################################
# Quantify transcripts using Salmon #
#####################################

# Create and enter working directory:
mkdir -p /home/workspace/ccorreia/globin/salmon_quant/human
cd !$

# Quantify transcripts from one FASTQ file to check if it works well:
salmon quant -i /home/workspace/ccorreia/globin/transcriptomes/human/human_index \
-l A -1 /home/workspace/ccorreia/globin/fastq_sequence/human/renamed_fastq/GD_TRUE_Subj12_L1_1.fastq.gz \
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
echo "salmon quant -i /home/workspace/ccorreia/globin/transcriptomes/human/human_index \
-l A -1 $file $file2 -2 $read1 $read2 \
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
echo "salmon quant -i /home/workspace/ccorreia/globin/transcriptomes/human/human_index \
-l A -1 $file $file2 -2 $read1 $read2 \
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
for file in `find /home/workspace/ccorreia/globin/salmon_quant/human \
-name quant.sf`; \
do oldname=`basename $file`; \
newname=`dirname $file | perl -p -e 's/.+(GD|NGD_.+_Subj\d+)/$1/'`; \
path=`dirname $file`; \
mv $file $path/${newname}_$oldname; \
done

# Move all *quant.sf files to a temporary folder:
mkdir /home/workspace/ccorreia/globin/salmon_quant/human/human_TPM
for file in `find /home/workspace/ccorreia/globin/salmon_quant/human \
-name *_Subj*quant.sf`; \
do cp $file -t /home/workspace/ccorreia/globin/salmon_quant/human/human_TPM; \
done

# Transfer all files from Rodeo to laptop:
scp -r ccorreia@servername:/home/workspace/ccorreia/globin/salmon_quant/human/human_TPM .

# Remove tmp folder from Rodeo:
rm -r human_TPM

# Append sample name to all log files:
for file in `find /home/workspace/ccorreia/globin/salmon_quant/human/ \
-name salmon_quant.log`; \
do oldname=`basename $file`; \
newname=`dirname $file | perl -p -e 's/.+(.+D_.+_Subj\d+).+/$1/'`; \
path=`dirname $file`; \
mv $file $path/${newname}_$oldname; \
done

# Gather salmon log information from all samples into one file:
for file in `find /home/workspace/ccorreia/globin/salmon_quant/human/ \
-name *D_*salmon_quant.log`; \
do echo echo \
"\`basename $file\` \
\`grep 'likely library type' $file | awk '{print \$12}'\` \
\`grep 'total fragments' $file | awk '{print \$2}'\` \
\`grep 'total reads' $file | awk '{print \$6}'\` \
\`grep 'Mapping rate' $file | awk '{print \$8}'\` >> \
salmon_summary_human.txt" >> salmon_summary_human.sh
done

chmod 755 salmon_summary_human.sh
./salmon_summary_human.sh

sed -i $'1 i\\\nFile_name Library_type Total_fragments Total_reads Mapping_rate(%)' \
salmon_summary_human.txt

# Transfer summary file from Rodeo to laptop:
scp -r ccorreia@remoteserver:/home/workspace/ccorreia/globin/salmon_quant/human/salmon_summary_human.txt .

# Following steps will be performed in R.


