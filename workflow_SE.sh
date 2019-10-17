#!/bin/bash
# -------------------------------------------------------

#Start of Pipline workflow_SE

#: <<'END'
echo [`date +"%Y-%m-%d %H:%M:%S"`]"----------Creating Directory hierachy------"
mkdir -p ~/data/workflow_SE/results/fastqc/
mkdir -p ~/data/workflow_SE/results/fastq_files/
mkdir -p ~/data/workflow_SE/results/cutadapt/
mkdir -p ~/data/workflow_SE/results/featureCounts/
mkdir -p ~/data/workflow_SE/results/hisat2/
mkdir -p ~/data/workflow_SE/results/samtools/

echo [`date +"%Y-%m-%d %H:%M:%S"`]"----------Made the Directories-------"
#END

#SRA_toolkit- Fastq--dump
#: << 'END'
echo [`date +"%Y-%m-%d %H:%M:%S"`] "----STARTING THE PIPELINE------"
echo "---------Running fastq-dump-----------"

cd ~/data/workflow_SE/sra/

##SRA to Fastq files


for file in $(ls | sed 's/.sra//')
do

 fastq-dump -F --gzip $file.sra

done

mv *.fastq.gz ~/data/workflow_SE/results/fastq_files/

echo "-----Moved fastq to results fastq_files folder------"


#END




##------------------------------------------------##
#: << 'END'

## Quality Check using Fastqc##

echo [`date +"%Y-%m-%d %H:%M:%S"`]"----------Doing Trimming and Quality check------"

cd ~/data/workflow_SE/results/fastq_files

#find ~/data/workflow_SE/results/fastq_files -name "*.fastq.gz" | sort | paste - - | while read A B 

echo "Doing quality checking"

for file in $(ls)
do 
  fastqc ${file}
done

mv *fastqc.html ~/data/workflow_SE/results/fastqc/
mv *fastqc.zip ~/data/workflow_SE/results/fastqc/

echo [`date +"%Y-%m-%d %H:%M:%S"`] "------------Done quality check and report generated!!----------"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "------------Moved to results fastqc!!----------"
#END

#--------------------------------------------###

#: << 'END'

##Trimming by Cutadapt
echo [`date +"%Y-%m-%d %H:%M:%S"`] "------Trimming by cutadapt------"

cd ~/data/workflow_SE/results/fastq_files/

for file in $(ls | sed 's/.fastq.gz//')

do

  cutadapt -m 10 -q 20 -j 8 -o ${file}\_trimmed.fastq.gz ${file}.fastq.gz

#USE your required adapter after -a

done
 
echo [`date +"%Y-%m-%d %H:%M:%S"`] "------Done trimming-----"


mv *_trimmed.fastq.gz ~/data/workflow_SE/results/cutadapt/

echo [`date +"%Y-%m-%d %H:%M:%S"`] "---------Moved to results cutadapt--------"

#END

#: << 'END'

##Index building and Read alignment using hisat2## 

cd ~/data/workflow_SE/reference_genome/

echo [`date +"%Y-%m-%d %H:%M:%S"`]"-----Building indices-----"
genome=~/data/workflow_SE/reference_genome/all.con.fa

##building index
hisat2-build $genome index 
#END

##----------------------------------------##

#: << 'END'

##Doing Alignment
echo [`date +"%Y-%m-%d %H:%M:%S"`]"----------Aligning with indices--------"
cd ~/data/workflow_SE/results/cutadapt/

for file in $(ls)
 
 do

  hisat2 -p 8 --dta -x ~/data/workflow_SE/reference_genome/index -U ${file} -S  ~/data/workflow_SE/results/hisat2/${file}.sam 
 

 done


 echo [`date +"%Y-%m-%d %H:%M:%S"`]"--------Done alignment and moved SAM files in results hisat2------"


#END

##-----------------------------------------##

#: << 'END'

##Converting sam files to bam files using SAMtools

echo [`date +"%Y-%m-%d %H:%M:%S"`]"------------Running SAM Tools---------"

cd ~/data/workflow_SE/results/hisat2/

for file in $(ls)
do
 samtools sort -@ 8 -o ${file}.bam ${file}
done

echo "-----converted sam to bam-----"

mv *.bam ~/data/workflow_SE/results/samtools/

echo [`date +"%Y-%m-%d %H:%M:%S"`]"--------Moved BAM file to samtool folder of results-----------"

echo [`date +"%Y-%m-%d %H:%M:%S"`]"--------Done with SAM tools---------"
#END

## ---------------------------------------------------------##

#: << 'END'

###FeatureCount tool
#Path to gtf file
gtf=~/data/workflow_SE/gtf/all.gtf

echo [`date +"%Y-%m-%d %H:%M:%S"`]"---------------Generating feature counts-------------"
 cd ~/data/workflow_SE/results/samtools/
  for file in $(ls)

 do

 featureCounts -t exon -g ID -a $gtf -o counts2.txt -M *.bam

 done

echo "---------Done Generating count data-----------"

mv *.txt ~/data/workflow_SE/results/featureCounts/
mv *.summary ~/data/workflow_SE/results/featureCounts/

echo "--------Results are in featureCount folder of results---------"
echo [`date +"%Y-%m-%d %H:%M:%S"`]"--------Finished with featurecounts-----"

#END

##-----------------------------------------------##
echo [`date +"%Y-%m-%d %H:%M:%S"`] "--------END OF PIPELINE---------"
exit 
