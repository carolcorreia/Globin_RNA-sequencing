##############################################################################
#      RNA-seq analysis of reticulocyte-derived globin gene transcripts      #
#             (HBA and HBB) in bovine peripheral blood samples               #
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
mkdir -p /home/workspace/ccorreia/globin/transcriptomes/cattle/source_file
cd !$

# Download and unzip the transcriptome:
nohup wget ftp://ftp.ensembl.org/pub/release-88/fasta/bos_taurus/cdna/Bos_taurus.UMD3.1.cdna.all.fa.gz &
gunzip Bos_taurus.UMD3.1.cdna.all.fa.gz


##############################################
# Build the transcriptome index using Salmon #
##############################################

# Required software is Salmon 0.8.2, consult manual/tutorial for details:
http://salmon.readthedocs.io/en/latest/

# Enter working directory:
cd /home/workspace/ccorreia/globin/transcriptomes/cattle

# Build an index for quasi-mapping:
nohup salmon index -t \
/home/workspace/ccorreia/globin/transcriptomes/cattle/source_file/Bos_taurus.UMD3.1.cdna.all.fa \
-i cattle_index --type quasi -k 31 -p 20 &


#####################################
# Quantify transcripts using Salmon #
#####################################

# Create and enter working directory:
mkdir -p /home/workspace/ccorreia/globin/salmon_quant/cattle
cd !$

# Quantify transcripts from one FASTQ file to check if it works well:
salmon quant -i /home/workspace/ccorreia/globin/transcriptomes/cattle/cattle_index \
-l A -1 \
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
/home/workspace/ccorreia/globin/transcriptomes/cattle/cattle_index \
-l A -1 $file $file2 $file3 $file4 $file5 \
-2 $read1 $read2 $read3 $read4 $read5 \
-p 15 -o ./$sample" \
>> quantify.sh; \
done

# Split and run all scripts on Rodeo:
split -d -l 5 quantify.sh quantify.sh.
for script in `ls quantify.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Append sample name to all quant.sf files:
for file in `find /home/workspace/ccorreia/globin/salmon_quant/cattle \
-name quant.sf`; \
do oldname=`basename $file`; \
newname=`dirname $file | perl -p -e 's/.+(A\d\d\d\d_W\-1_F)/$1/'`; \
path=`dirname $file`; \
mv $file $path/${newname}_$oldname; \
done

# Move all *quant.sf files to a temporary folder:
mkdir /home/workspace/ccorreia/globin/salmon_quant/cattle/cattle_TPM
for file in `find /home/workspace/ccorreia/globin/salmon_quant/cattle \
-name A*quant.sf`; \
do cp $file -t /home/workspace/ccorreia/globin/salmon_quant/cattle/cattle_TPM; \
done

# Transfer all files from Rodeo to laptop:
scp -r ccorreia@servername:/home/workspace/ccorreia/globin/salmon_quant/cattle/cattle_TPM .

# Remove tmp folder from Rodeo:
rm -r cattle_TPM

# Append sample name to all log files:
for file in `find /home/workspace/ccorreia/globin/salmon_quant/cattle/ \
-name salmon_quant.log`; \
do oldname=`basename $file`; \
newname=`dirname $file | perl -p -e 's/.+(A\d\d\d\d_W\-1_F).+/$1/'`; \
path=`dirname $file`; \
mv $file $path/${newname}_$oldname; \
done

# Gather salmon log information from all samples into one file:
for file in `find /home/workspace/ccorreia/globin/salmon_quant/cattle/ \
-name 7*salmon_quant.log`; \
do echo echo \
"\`basename $file\` \
\`grep 'likely library type' $file | awk '{print \$12}'\` \
\`grep 'total fragments' $file | awk '{print \$2}'\` \
\`grep 'total reads' $file | awk '{print \$6}'\` \
\`grep 'Mapping rate' $file | awk '{print \$8}'\` >> \
salmon_summary_cattle.txt" >> salmon_summary_cattle.sh
done

chmod 755 salmon_summary_cattle.sh
./salmon_summary_cattle.sh

sed -i $'1 i\\\nFile_name Library_type Total_fragments Total_reads Mapping_rate(%)' \
salmon_summary_cattle.txt

# Transfer summary file from Rodeo to laptop:
scp -r ccorreia@servername:/home/workspace/ccorreia/globin/salmon_quant/cattle/salmon_summary_cattle.txt .

# Following steps were performed in R.










