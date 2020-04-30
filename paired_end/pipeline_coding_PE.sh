#################################################################
#################     USER SCRIPT START     #####################
#################################################################

start_time=`date +%s`
echo "start time:" $start_time

#################################################################
#################   DECLARATIONS / PATHS   ######################
#################################################################
proj="rice_rsolani_zheng1" #manually set job name here before running, one project can have multiole jobs
job="rice_rsolani_zheng1" #manually set job name here before running
sra="/home/cluster/ankur/project/rice_rsolani_zheng/data/sra" #sra files (source) directory (no trailing /)
fastq_f="" #fastq files (source) directory (no trailing /)
ref_g="/home/cluster/ankur/project/rice_rsolani_zheng/data/ref/ref_sativa.fa" #full path of your (source) reference genome (.fa)
gtf_f="/home/cluster/ankur/project/rice_rsolani_zheng/data/ref/Oryza_sativa.IRGSP.gtf" #full path of (source) gtf (.gtf)
deseq2="/path/to/deseq2/R/script.R" #full path of the DeSeq2 R script (.R)

#################################################################
################ Start of Pipline workflow_PE ###################
#################################################################

#: <<'END'
echo [`date +"%Y-%m-%d %H:%M:%S"`]"----------Creating Directory hierachy------"
mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastqc/
mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastq_files/
mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastp/
mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/featureCounts/
mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/hisat2/
mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/samtools/

mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/sra/	#chiranjeevdas
mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/reference_genome/	#chiranjeevdas
mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/gtf/	#chiranjeevdas
mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastp_preqc/	#chiranjeevdas
mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastp_reports/
mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastqc_post_fastp/	#chiranjeevdas

#mkdir -p ~/projects/$proj/pipeline_result/$job/scripts/ #chiranjeevdas

echo [`date +"%Y-%m-%d %H:%M:%S"`]"----------Made the Directories-------"
#END

#################################################################
############## ADDITIONAL TASKS/ BYPASSES #chiranjeevdas#########
#################################################################

#KEEP THIS AREA CLEAR WHEN NOT IN USE

#: <<'END'

echo [`date +"%Y-%m-%d %H:%M:%S"`]"----------Doing additional tasks-------"
rsync -a -h -v -r -P -t $fastq_f/*.gz ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastq_files/ 	#chiranjeevdas
echo [`date +"%Y-%m-%d %H:%M:%S"`]"----------Additional tasks complete-------"

#END

#################################################################

#SRA_toolkit- Fastq--dump

: << 'END'
echo [`date +"%Y-%m-%d %H:%M:%S"`] "----STARTING THE PIPELINE------"
echo "---------Running fastq-dump-----------"

rsync -a -h -v -r -P -t $sra/* ~/projects/$proj/pipeline_result/$job/data/workflow_PE/sra/	#chiranjeevdas
cd ~/projects/$proj/pipeline_result/$job/data/workflow_PE/sra/

##SRA to Fastq files

for file in $(ls | sed 's/.sra//')
do

 fasterq-dump $file.sra --split-files -v -p -b 100 -c 1024 -m 10240 -e 50 	#resource limit #chiranjeevdas

done

pigz -v *.fastq #multhreaded gzip #resource limit #chiranjeevdas

rsync -a -h -v -r -P -t --remove-source-files *.fastq.gz ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastq_files/


echo "-----Moved fastq to results fastq_files folder------"


END


##------------------------------------------------##
: << 'END'

## Quality Check using Fastqc##

echo [`date +"%Y-%m-%d %H:%M:%S"`]"----------Doing Quality check------"

cd ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastq_files

echo "Doing quality checking"

f_ls="$(ls)"

fastqc -t 50 $f_ls #resource limit #chiranjeevdas
#-t = number of threads (250 MB memory required per thread)

rsync -a -h -v -r -P -t --remove-source-files *fastqc.html ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastqc/
rsync -a -h -v -r -P -t --remove-source-files *fastqc.zip ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastqc/

echo [`date +"%Y-%m-%d %H:%M:%S"`] "------------Done quality check and report generated!!----------"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "------------Moved to results fastqc!!----------"
END

#--------------------------------------------###

#: << 'END'

##Trimming by fastp
echo [`date +"%Y-%m-%d %H:%M:%S"`] "------Trimming by fastp------"
    
cd ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastq_files/

find ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastq_files/ -name "*.fastq.gz" | sort | paste - - | while read A B
do

fastp -V --thread=16 --length_required=10 --qualified_quality_phred=30 --in1=${A}.fastq.gz --in2=${B}.fastq.gz --out1=${A}\_trimmed.fastq.gz --out2=${A}\_trimmed.fastq.gz --json=${A}.json --html=${A}.html	#chiranjeevdas

rsync -a -h -v -r -P -t --remove-source-files *.html *.json ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastp_reports/
  #--thread= number of worker threads (max 16)
  #USE your required adapter after -a, default = automatic detection

done

echo [`date +"%Y-%m-%d %H:%M:%S"`] "------Done trimming-----"

rsync -a -h -v -r -P -t  *_trimmed.fastq.gz ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastp_preqc/
rsync -a -h -v -r -P -t --remove-source-files *_trimmed.fastq.gz ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastp/

echo [`date +"%Y-%m-%d %H:%M:%S"`] "---------Moved to results fastp--------"

#END

#: << 'END'

## Quality Check using Fastqc AGAIN##

echo [`date +"%Y-%m-%d %H:%M:%S"`]"----------Doing Quality check after trimming------"

cd ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastp_preqc/

#find ~/data/workflow_PE/results/fastq_files -name "*.fastq.gz" | sort | paste - - | while read A B 

echo "Doing quality checking"

f_ls="$(ls)"

fastqc -t 50 $f_ls #resource limit #chiranjeevdas
#-t = number of threads (250 MB memory required per thread)

rsync -a -h -v -r -P -t --remove-source-files *fastqc.html ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastqc_post_fastp/
rsync -a -h -v -r -P -t --remove-source-files *fastqc.zip ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastqc_post_fastp/

echo [`date +"%Y-%m-%d %H:%M:%S"`] "------------Done quality check and report generated!!----------"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "------------Moved to results fastqc_post_fastp!!----------"


#END


#: << 'END'

##Index building and Read alignment using hisat2## 

rsync -a -h -v -r -P -t  $ref_g ~/projects/$proj/pipeline_result/$job/data/workflow_PE/reference_genome/ref_sativa.fa #chiranjeevdas
cd ~/projects/$proj/pipeline_result/$job/data/workflow_PE/reference_genome/ #chiranjeevdas

echo [`date +"%Y-%m-%d %H:%M:%S"`]"-----Building indices-----"
genome=~/projects/$proj/pipeline_result/$job/data/workflow_PE/reference_genome/ref_sativa.fa

##building index
hisat2-build -p 50 $genome index #resource limit #chiranjeevdas 
#-p = number of threads

#END

##----------------------------------------##

#: << 'END'

##Doing Alignment
echo [`date +"%Y-%m-%d %H:%M:%S"`]"----------Aligning with indices--------"
cd ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastp/

find ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastp/ -name "*.gz_trimmed_fastq.gz" | sort | paste - - | while read A B
 
 do

hisat2 --threads 50 --dta -x ~/projects/$proj/pipeline_result/$job/data/workflow_PE/reference_genome/index -1 ${A} -2 ${B} -S  ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/hisat2/${A}.sam		#resource limit #chiranjeevdas 
 	#--threads = number of simultaneous alignments

 done

echo [`date +"%Y-%m-%d %H:%M:%S"`]"--------Done alignment and moved SAM files in results hisat2------"

#END

##-----------------------------------------##

#: << 'END'

##Converting sam files to bam files using SAMtools

echo [`date +"%Y-%m-%d %H:%M:%S"`]"------------Running SAM Tools---------"

cd ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/hisat2/

for file in $(ls)
do
 samtools sort -@ 30 -o ${file}.bam ${file} 	#resource limit #chiranjeevdas
	#-@ number of threads (in addition to main thread)
done

echo "-----converted sam to bam-----"

rsync -a -h -v -r -P -t --remove-source-files *.bam ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/samtools/

echo [`date +"%Y-%m-%d %H:%M:%S"`]"--------Moved BAM file to samtool folder of results-----------"

echo [`date +"%Y-%m-%d %H:%M:%S"`]"--------Done with SAM tools---------"

#END

## ---------------------------------------------------------##


###FeatureCount tool

echo [`date +"%Y-%m-%d %H:%M:%S"`]"---------------Generating feature counts-------------"

cd ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/samtools/

featureCounts -T 50 -t gene -g gene_id -a $gtf_f -o counts.txt -M *.bam	#resource limit #chiranjeevdas
	#-T number of threads

echo "---------Done Generating count data-----------"

rsync -a -h -v -r -P -t --remove-source-files *.txt ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/featureCounts/
rsync -a -h -v -r -P -t --remove-source-files *.summary ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/featureCounts/

echo "--------Results are in featureCount folder of results---------"

echo [`date +"%Y-%m-%d %H:%M:%S"`]"--------Finished with featurecounts-----"

#END

##-----------------------------------------##

: << 'END'

### DeSeq2 in R

echo [`date +"%Y-%m-%d %H:%M:%S"`]"---------------Calling R for DeSeq2-------------"

cp $deseq2 ~/projects/$proj/pipeline_result/$job/scripts/deseq2.R
cd ~/projects/$proj/pipeline_result/$job/scripts/

r deseq2.R

mv res.csv bak_res.csv
echo -n "", > res.csv; cat bak_res.csv >> res.csv #fixes the left shift of column names
rm bak_res.csv

mv res_PAdj_cutoff.csv bak_res_PAdj_cutoff.csv
echo -n "", > res_PAdj_cutoff.csv; cat bak_res_PAdj_cutoff.csv >> res_PAdj_cutoff.csv #fixes the left shift of column names
rm bak_res_PAdj_cutoff.csv 

echo [`date +"%Y-%m-%d %H:%M:%S"`]"---------------DeSeq2 Complete-------------"

END

##-----------------------------------------------##
echo [`date +"%Y-%m-%d %H:%M:%S"`] "--------END OF PIPELINE---------"
############################################################

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

############################################################

exit

############################################################
#################### END OF SCRIPT #########################
############################################################
