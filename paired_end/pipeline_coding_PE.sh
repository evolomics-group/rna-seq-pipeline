#!/bin/bash
#SBATCH --job-name=45
#SBATCH --ntasks=50

#################################################################
#################     USER INSTRUCTIONS     #####################
#################################################################

# 1. Modify the "job-name" above (OPTIONAL)
# 2. In the "DECLARATIONS/PATHS" section below, put a project name.
# 3. In the "DECLARATIONS/PATHS" section below, put a unique job name. NB. This is specific to each run of the pipeline and must be edited every time you run it.
# 4. Populate the 'sra', 'fatstq_f', 'ref_g', 'gtf_f' & 'deseq2' paths (you may ommit paths if you are not running that step). Ensure the paths/files exist.
# 5. Check the number of threads you CPU can process and change the variables in the lines which can be found by searching for the tag "#resource limit"
# 5. Save the edited script as a copy and run with sbatch <script name>

#################################################################
#################     USER SCRIPT START     #####################
#################################################################

start_time=`date +%s`
echo "start time:" $start_time

#################################################################
#################   DECLARATIONS / PATHS   ######################
#################################################################
proj="rsolani" #manually set job name here before running, one project can have multiole jobs
job="project_y" #manually set job name here before running

sra="/home/cluster/akash/projects/rsolani_nopaper/data/srarawdata" #sra files (source) directory (no trailing /)
fastq_f="path/to/fastq/folder" #fastq files (source) directory (no trailing /)
ref_g="/home/cluster/evolomics/projects/rice_rsolani/data/ref_genome/ref_chr_mt_pt_o_sativa.fa" #full path of your (source) reference genome file ending in .fa OR .fasta
gtf_f="/home/cluster/evolomics/projects/rice_rsolani/data/gtf/Oryza_sativa.IRGSP-1.0.40.gtf" #full path of (source) gtf file ending in .gtf
deseq2="/path/to/deseq2/R/script.R" #full path of the DeSeq2 R script (.R)

#################################################################
################ Start of Pipline workflow_PE ###################
#################################################################
echo [`date +"%Y-%m-%d %H:%M:%S"`] "----STARTING THE PIPELINE------"

#: <<'END'
echo [`date +"%Y-%m-%d %H:%M:%S"`] "----------Creating Directory hierachy------"
mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastqc/
mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastq_files/
mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastp/
mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastp_reports/
mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/featureCounts/
mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/hisat2/
mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/samtools/

mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/reference_genome/
mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/gtf/

#: << 'END'

#ADD A HASH BEFORE THE COLON ABOVE IF INPUT FILES ARE SRA (ALSO FIND THE CORRESPONDING TAG BELOW)
#REMOVE THE HASH (IF ANY) BEFORE THE COLON ABOVE IF INPUT FILES ARE FASTQ/FASTQ.GZ (ALSO FIND THE CORRESPONDING TAG BELOW)

echo [`date +"%Y-%m-%d %H:%M:%S"`] "---------Running fastq-dump-----------"

cd ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastq_files/

##SRA to Fastq files

for file in $(ls $sra)
do

fasterq-dump $sra/$file --split-files -v -p -b 100 -c 1024 -m 10240 -e 50       #resource limit

done

echo [`date +"%Y-%m-%d %H:%M:%S"`] "---------Conversion to FastQ complete. Compressing now.-----------"

pigz -v *.fastq #multhreaded gzip #resource limit

fastq_f=~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastq_files/

echo [`date +"%Y-%m-%d %H:%M:%S"`] "-----All SRAs split and converted to fastq.gz----"

#END

#ADD A HASH BEFORE THE END ABOVE IF INPUT FILES ARE SRA (ALSO FIND THE CORRESPONDING TAG ABOVE)
#REMOVE THE HASH (IF ANY) BEFORE THE END ABOVE IF INPUT FILES ARE FASTQ/FASTQ.GZ (ALSO FIND THE CORRESPONDING TAG ABOVE)


##Trimming by fastp
echo [`date +"%Y-%m-%d %H:%M:%S"`] "------Trimming by fastp------"

cd ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastp
find $fastq_f -name "*.fastq.gz" | sort | paste - - | while read A B

do

a=`basename ${A} | sed 's/.sra_1/_1/' | awk -F "." '{print $1}'`
b=`basename ${B} | sed 's/.sra_2/_2/' | awk -F "." '{print $1}'`

echo ""
echo "Processing $a and $b"

fastp --thread=16 --length_required=10 --qualified_quality_phred=32 --in1=${A} --in2=${B} --out1=$a\_trimmed.fastq.gz --out2=$b\_trimmed.fastq.gz --json=$a.json --html=$a.html
  #--thread= number of worker threads (max 16)
  #USE your required adapter after -a, default = automatic detection

mv -v $a.json ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastp_reports/`echo $a | awk -F "_" '{print $1".json"}'`
mv -v $a.html ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastp_reports/`echo $a | awk -F "_" '{print $1".html"}'`

echo ""
echo "$a and $b trimmed. Report generated and moved to results fastp_reports!!!"

done

echo [`date +"%Y-%m-%d %H:%M:%S"`] "------All trimming completed!!!----"

echo [`date +"%Y-%m-%d %H:%M:%S"`]"----------Doing QC after trimming------"

cd ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastp/

f_ls="$(ls)"

fastqc -q -t 50 $f_ls #resource limit
#-t = number of threads (250 MB memory required per thread)

ls | awk -F "." '{print $1}' | uniq | while read report

do

mv -v $report*.html ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastqc/`echo $report | awk -F "." '{print $1}'`.html
mv -v $report*.zip ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastqc/`echo $report | awk -F "." '{print $1}'`.zip

done

echo [`date +"%Y-%m-%d %H:%M:%S"`] "------------Done quality check and report generated!!----------"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "------------Moved to results fastqc!!----------"


#END

#: << 'END'

echo [`date +"%Y-%m-%d %H:%M:%S"`] "------------Preparing for Alignment----------"

##Index building and Read alignment using hisat2## 

echo [`date +"%Y-%m-%d %H:%M:%S"`] "------------Locating and placing reference genome------------"
cp -v $ref_g ~/projects/$proj/pipeline_result/$job/data/workflow_PE/reference_genome/ref_genome.fa
cd ~/projects/$proj/pipeline_result/$job/data/workflow_PE/reference_genome/

echo [`date +"%Y-%m-%d %H:%M:%S"`] "-----Building indices-----"

##building index
hisat2-build -p 50 ~/projects/$proj/pipeline_result/$job/data/workflow_PE/reference_genome/ref_genome.fa index #resource limit
#-p = number of threads

#END

##Doing Alignment
echo [`date +"%Y-%m-%d %H:%M:%S"`] "----------Aligning with indices--------"

cd ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastp/

find ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastp/ -name "*_trimmed.fastq.gz" | sort | paste - - | while read A B

do

a=`basename ${A} | awk -F "." '{print $1}' | awk -F "_" '{print $1}'`

hisat2 --threads 50 --dta -x ~/projects/$proj/pipeline_result/$job/data/workflow_PE/reference_genome/index -1 ${A} -2 ${B} -S  ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/hisat2/$a.sam		#resource limit #chiranjeevdas 
 	#--threads = number of simultaneous alignments
done

echo [`date +"%Y-%m-%d %H:%M:%S"`]"--------Done alignment and placed SAM files in results hisat2------"

#END

##-----------------------------------------##

#: << 'END'

##Converting sam files to bam files using SAMtools

echo [`date +"%Y-%m-%d %H:%M:%S"`] "------------Running SAM Tools to Convert SAM to BAM---------"

cd ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/hisat2/

for file in $(ls)

do

a=`echo ${file} | awk -F "." '{print $1}'`

echo "Processing ${file}"
samtools sort -@ 30 -o $a.bam ${file} 	#resource limit
	#-@ number of threads (in addition to main thread)
echo "${file} converted"
done

echo [`date +"%Y-%m-%d %H:%M:%S"`] "-----Converted all SAM files to BAM-----"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "--------Moving BAM files to results samtools-----------"

mv -v *.bam ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/samtools/

echo [`date +"%Y-%m-%d %H:%M:%S"`] "--------Moved BAM files-----------"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "--------Done with SAM tools !!!---------"

#END

## ---------------------------------------------------------##


###FeatureCount tool

echo [`date +"%Y-%m-%d %H:%M:%S"`] "---------------Generating feature counts...-------------"

cd ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/samtools/

featureCounts -T 50 -t gene -g gene_id -a $gtf_f -o counts.txt -M *.bam	#resource limit
	#-T number of threads

echo [`date +"%Y-%m-%d %H:%M:%S"`] "---------Done Generating count data-----------"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "---------Moving Results...-----------"

mv -v *.txt ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/featureCounts/
mv -v *.summary ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/featureCounts/

echo [`date +"%Y-%m-%d %H:%M:%S"`] "--------Results moved to results featureCounts---------"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "--------Finished with featureCounts!!!-----"

#END

##-----------------------------------------##

: << 'END'
### DeSeq2 in R

echo [`date +"%Y-%m-%d %H:%M:%S"`] "---------------Calling R for DeSeq2-------------"

cp $deseq2 ~/projects/$proj/pipeline_result/$job/scripts/deseq2.R
cd ~/projects/$proj/pipeline_result/$job/scripts/

r deseq2.R

mv -v res.csv bak_res.csv
echo -n "", > res.csv; cat bak_res.csv >> res.csv #fixes the left shift of column names
rm bak_res.csv

mv -v res_PAdj_cutoff.csv bak_res_PAdj_cutoff.csv
echo -n "", > res_PAdj_cutoff.csv; cat bak_res_PAdj_cutoff.csv >> res_PAdj_cutoff.csv #fixes the left shift of column names
rm bak_res_PAdj_cutoff.csv

echo [`date +"%Y-%m-%d %H:%M:%S"`]"---------------DeSeq2 Complete-------------"

END

##-----------------------------------------------##

echo [`date +"%Y-%m-%d %H:%M:%S"`]"---------------Cleaning Up-------------"

echo [`date +"%Y-%m-%d %H:%M:%S"`]"---------------Deleting SAM files-------------"
rm -v ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/hisat2/*
echo [`date +"%Y-%m-%d %H:%M:%S"`]"---------------Deleting Trimmed FastQ Files-------------"
rm -v ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastp/*
echo [`date +"%Y-%m-%d %H:%M:%S"`]"---------------Deleting FastQ files-------------"
rm -v ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastq_files/*
echo [`date +"%Y-%m-%d %H:%M:%S"`]"---------------Clean Up Completed !!!-------------"

echo [`date +"%Y-%m-%d %H:%M:%S"`]"---------------Compressing BAM files-------------"
cd ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/samtools/
echo -e "\nCompressing:\n"
pigz -v -9 *
echo [`date +"%Y-%m-%d %H:%M:%S"`]"---------------Compression Completed !!!-------------"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "--------END OF PIPELINE---------"

end_time=`date +%s`
echo "############################################################"
echo
echo "time taken:"
echo
echo "start time:" $start_time
echo "end time:" $end_time
echo
run_time=$((end_time-start_time))
echo "runtime (in seconds):" $run_time
echo -n "runtime (in minutes): "; awk "BEGIN {print $run_time/60}"
echo -n "runtime (in hours): "; awk "BEGIN {print $run_time/3600}"
echo
echo "############################################################"

exit

############################################################
#################### END OF SCRIPT #########################
############################################################
#edited_on_github_atom
