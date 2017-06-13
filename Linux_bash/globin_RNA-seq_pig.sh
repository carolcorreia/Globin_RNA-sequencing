##############################################################################
#      RNA-seq analysis of reticulocyte-derived globin gene transcripts      #
#            (HBA and HBB) in porcine peripheral blood samples               #
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
mkdir -p /home/workspace/ccorreia/globin/transcriptomes/pig/source_file
cd !$

# Download and unzip the transcriptome:
nohup wget ftp://ftp.ensembl.org/pub/release-88/fasta/sus_scrofa/cdna/Sus_scrofa.Sscrofa10.2.cdna.all.fa.gz &
gunzip Sus_scrofa.Sscrofa10.2.cdna.all.fa.gz


##############################################
# Build the transcriptome index using Salmon #
##############################################

# Required software is Salmon 0.8.2, consult manual/tutorial for details:
http://salmon.readthedocs.io/en/latest/

# Enter working directory:
cd /home/workspace/ccorreia/globin/transcriptomes/pig

# Build an index for quasi-mapping:
nohup salmon index -t \
/home/workspace/ccorreia/globin/transcriptomes/pig/source_file/Sus_scrofa.Sscrofa10.2.cdna.all.fa \
-i pig_index --type quasi -k 31 -p 20 &


#####################################
# Quantify transcripts using Salmon #
#####################################

# Create and enter working directory:
mkdir -p /home/workspace/ccorreia/globin/salmon_quant/pig
cd !$

# Quantify transcripts from one FASTQ file to check if it works well:
salmon quant -i /home/workspace/ccorreia/globin/transcriptomes/pig/pig_index \
-l A -1 /home/workspace/ccorreia/globin/fastq_sequence/pig/HI.0751.004.Index_12.GCswine-5037-28DPI-WB-7413C/trimmed_HI.0751.004.Index_12.GCswine-5037-28DPI-WB-7413C-mRNA_R1.fastq.gz \
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
echo "salmon quant -i /home/workspace/ccorreia/globin/transcriptomes/pig/pig_index \
-l A -1 $file -2 $file2 \
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
echo "salmon quant -i /home/workspace/ccorreia/globin/transcriptomes/pig/pig_index \
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
for file in `find /home/workspace/ccorreia/globin/salmon_quant/pig \
-name quant.sf`; \
do oldname=`basename $file`; \
newname=`dirname $file | perl -p -e 's/.+(\d\d\d\d.+)/$1/'`; \
path=`dirname $file`; \
mv $file $path/${newname}_$oldname; \
done

# Move all *quant.sf files to a temporary folder:
mkdir /home/workspace/ccorreia/globin/salmon_quant/pig/pig_TPM
for file in `find /home/workspace/ccorreia/globin/salmon_quant/pig \
-name 7*quant.sf`; \
do cp $file -t /home/workspace/ccorreia/globin/salmon_quant/pig/pig_TPM; \
done

# Transfer all files from Rodeo to laptop:
scp -r ccorreia@servername:/home/workspace/ccorreia/globin/salmon_quant/pig/pig_TPM .

# Remove tmp folder from Rodeo:
rm -r pig_TPM

# Append sample name to all log files:
for file in `find /home/workspace/ccorreia/globin/salmon_quant/pig/ \
-name salmon_quant.log`; \
do oldname=`basename $file`; \
newname=`dirname $file | perl -p -e 's/.+(7\d\d\d\w+).+/$1/'`; \
path=`dirname $file`; \
mv $file $path/${newname}_$oldname; \
done

# Gather salmon log information from all samples into one file:
for file in `find /home/workspace/ccorreia/globin/salmon_quant/pig/ \
-name 7*salmon_quant.log`; \
do echo echo \
"\`basename $file\` \
\`grep 'likely library type' $file | awk '{print \$12}'\` \
\`grep 'total fragments' $file | awk '{print \$2}'\` \
\`grep 'total reads' $file | awk '{print \$6}'\` \
\`grep 'Mapping rate' $file | awk '{print \$8}'\` >> \
salmon_summary_pig.txt" >> salmon_summary_pig.sh
done

chmod 755 salmon_summary_pig.sh
./salmon_summary_pig.sh

sed -i $'1 i\\\nFile_name Library_type Total_fragments Total_reads Mapping_rate(%)' \
salmon_summary_pig.txt

# Transfer summary file from Rodeo to laptop:
scp -r ccorreia@rodeo.ucd.ie:/home/workspace/ccorreia/globin/salmon_quant/pig/salmon_summary_pig.txt .

# Following steps will be performed in R.


