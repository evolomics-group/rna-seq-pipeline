#################################################################
#################     USER SCRIPT START     #####################
#################################################################

start_time=`date +%s` 
echo "start time:" $start_time

#################################################################
#################   DECLARATIONS / PATHS   ######################
#################################################################

job="rizoctonia8" #manually set job name here before running
#sra="/full/path/to/sra/files" #sra files (source) directory (no trailing /)
ref_g="/home/cluster/chiranjeev/sources/ref_chr_mt_pt_o_sativa.fa" #full path of your (source) reference genome (.fa)
gtf_f="/home/cluster/chiranjeev/sources/Oryza_sativa.IRGSP-1.0.40.gtf" #full path of (source) gtf (.gtf)
deseq2="/path/to/deseq2/R/script.R" #full path of the DeSeq2 R script (.R)

echo "job name:" $job
echo "sra location:" $sra
echo "reference genome location:" $ref_g
echo "gtf file location:" $gtf_f

#################################################################
################ Start of Pipline workflow_SE ###################
#################################################################

#: <<'END'
echo [`date +"%Y-%m-%d %H:%M:%S"`]"----------Creating Directory hierachy------"
mkdir -p ~/$job/data/workflow_SE/results/fastqc/
mkdir -p ~/$job/data/workflow_SE/results/fastq_files/
mkdir -p ~/$job/data/workflow_SE/results/fastp/
mkdir -p ~/$job/data/workflow_SE/results/featureCounts/
mkdir -p ~/$job/data/workflow_SE/results/hisat2/
mkdir -p ~/$job/data/workflow_SE/results/samtools/

mkdir -p ~/$job/data/workflow_SE/sra/	#chiranjeevdas
mkdir -p ~/$job/data/workflow_SE/reference_genome/	#chiranjeevdas
mkdir -p ~/$job/data/workflow_SE/gtf/	#chiranjeevdas
mkdir -p ~/$job/data/workflow_SE/results/fastp_preqc/	#chiranjeevdas
mkdir -p ~/$job/data/workflow_SE/results/fastp_reports/
mkdir -p ~/$job/data/workflow_SE/results/fastqc_post_fastp/	#chiranjeevdas

mkdir -p ~/$job/scripts/ #chiranjeevdas

mkdir -p ~/$job/data/workflow_SE/results/stringtie/gtfs/
mkdir -p ~/$job/data/workflow_SE/results/stringtie/abundance/
mkdir -p ~/$job/data/workflow_SE/results/stringtie/merge/gtfs/
mkdir -p ~/$job/data/workflow_SE/results/stringtie/merge/abundance/

echo [`date +"%Y-%m-%d %H:%M:%S"`]"----------Made the Directories-------"
#END

#################################################################
############## ADDITIONAL TASKS/ BYPASSES #chiranjeevdas#########
#################################################################

#KEEP THIS AREA CLEAR WHEN NOT IN USE

#: <<'END'

echo [`date +"%Y-%m-%d %H:%M:%S"`]"----------Doing additional tasks-------"
rsync -a -h -v -r -P -t /home/cluster/evolomics/common/data/rizoctonia/fastq/*.gz ~/$job/data/workflow_SE/results/fastq_files/ 	#chiranjeevdas
echo [`date +"%Y-%m-%d %H:%M:%S"`]"----------Additional tasks complete-------"

#END

#################################################################

#SRA_toolkit- Fastq--dump

: << 'END'
echo [`date +"%Y-%m-%d %H:%M:%S"`] "----STARTING THE PIPELINE------"
echo "---------Running fastq-dump-----------"
rsync -a -h -v -r -P -t $sra/* ~/$job/data/workflow_SE/sra/	#chiranjeevdas
cd ~/$job/data/workflow_SE/sra/
##SRA to Fastq files
for file in $(ls | sed 's/.sra//')
do
 fasterq-dump $file.sra -v -p -b 10 -c 100 -m 4096 -e 50 	#resource limit #chiranjeevdas
done
pigz -v *.fastq #multhreaded gzip #resource limit #chiranjeevdas
rsync -a -h -v -r -P -t --remove-source-files *.fastq.gz ~/data/workflow_SE/results/fastq_files/
echo "-----Moved fastq to results fastq_files folder------"
END


##------------------------------------------------##
#: << 'END'

## Quality Check using Fastqc##

echo [`date +"%Y-%m-%d %H:%M:%S"`]"----------Doing Quality check------"

cd ~/$job/data/workflow_SE/results/fastq_files

echo "Doing quality checking"

f_ls="$(ls)"

fastqc -t 50 $f_ls #resource limit #chiranjeevdas
#-t = number of threads (250 MB memory required per thread)

rsync -a -h -v -r -P -t --remove-source-files *fastqc.html ~/$job/data/workflow_SE/results/fastqc/
rsync -a -h -v -r -P -t --remove-source-files *fastqc.zip ~/$job/data/workflow_SE/results/fastqc/

echo [`date +"%Y-%m-%d %H:%M:%S"`] "------------Done quality check and report generated!!----------"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "------------Moved to results fastqc!!----------"
#END

#--------------------------------------------###

#: << 'END'

##Trimming by fastp
echo [`date +"%Y-%m-%d %H:%M:%S"`] "------Trimming by fastp------"

cd ~/$job/data/workflow_SE/results/fastq_files/

for file in $(ls | sed 's/.fastq.gz//')

do

  fastp -V --thread=16 --length_required=10 --qualified_quality_phred=30 --in1=${file}.fastq.gz --out1=${file}\_trimmed.fastq.gz --json=${file}.json --html=${file}.html	#chiranjeevdas
  rsync -a -h -v -r -P -t --remove-source-files *.html *.json ~/$job/data/workflow_SE/results/fastp_reports/
  #--thread= number of worker threads (max 16)
  #USE your required adapter after -a, default = automatic detection

done

echo [`date +"%Y-%m-%d %H:%M:%S"`] "------Done trimming-----"

rsync -a -h -v -r -P -t  *_trimmed.fastq.gz ~/$job/data/workflow_SE/results/fastp_preqc/
rsync -a -h -v -r -P -t --remove-source-files *_trimmed.fastq.gz ~/$job/data/workflow_SE/results/fastp/

echo [`date +"%Y-%m-%d %H:%M:%S"`] "---------Moved to results fastp--------"

#END

#: << 'END'

## Quality Check using Fastqc AGAIN##

echo [`date +"%Y-%m-%d %H:%M:%S"`]"----------Doing Quality check after trimming------"

cd ~/$job/data/workflow_SE/results/fastp_preqc/

#find ~/data/workflow_SE/results/fastq_files -name "*.fastq.gz" | sort | paste - - | while read A B 

echo "Doing quality checking"

f_ls="$(ls)"

fastqc -t 50 $f_ls #resource limit #chiranjeevdas
#-t = number of threads (250 MB memory required per thread)

rsync -a -h -v -r -P -t --remove-source-files *fastqc.html ~/$job/data/workflow_SE/results/fastqc_post_fastp/
rsync -a -h -v -r -P -t --remove-source-files *fastqc.zip ~/$job/data/workflow_SE/results/fastqc_post_fastp/

echo [`date +"%Y-%m-%d %H:%M:%S"`] "------------Done quality check and report generated!!----------"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "------------Moved to results fastqc_post_fastp!!----------"


#END


#: << 'END'

##Index building and Read alignment using hisat2## 

rsync -a -h -v -r -P -t  $ref_g ~/$job/data/workflow_SE/reference_genome/ref_genome.fa #chiranjeevdas
cd ~/$job/data/workflow_SE/reference_genome/ #chiranjeevdas

echo [`date +"%Y-%m-%d %H:%M:%S"`]"-----Building indices-----"
genome=~/$job/data/workflow_SE/reference_genome/ref_genome.fa

##building index
hisat2-build -p 50 $genome index #resource limit #chiranjeevdas 
#-p = number of threads

#END

##----------------------------------------##

#: << 'END'

##Doing Alignment
echo [`date +"%Y-%m-%d %H:%M:%S"`]"----------Aligning with indices--------"
cd ~/$job/data/workflow_SE/results/fastp/

for file in $(ls)
 
 do

  hisat2 --threads 50 --dta -x ~/$job/data/workflow_SE/reference_genome/index -U ${file} -S  ~/$job/data/workflow_SE/results/hisat2/${file}.sam		#resource limit #chiranjeevdas 
 	#--threads = number of simultaneous alignments

 done

echo [`date +"%Y-%m-%d %H:%M:%S"`]"--------Done alignment and moved SAM files in results hisat2------"

#END

##-----------------------------------------##

#: << 'END'

##Converting sam files to bam files using SAMtools

echo [`date +"%Y-%m-%d %H:%M:%S"`]"------------Running SAM Tools---------"

cd ~/$job/data/workflow_SE/results/hisat2/

for file in $(ls)
do
 samtools sort -@ 30 -o ${file}.bam ${file} 	#resource limit #chiranjeevdas
	#-@ number of threads (in addition to main thread)
done

echo "-----converted sam to bam-----"

rsync -a -h -v -r -P -t --remove-source-files *.bam ~/$job/data/workflow_SE/results/samtools/

echo [`date +"%Y-%m-%d %H:%M:%S"`]"--------Moved BAM file to samtool folder of results-----------"

echo [`date +"%Y-%m-%d %H:%M:%S"`]"--------Done with SAM tools---------"

#END

## ---------------------------------------------------------##

#: << 'END'

### String Tie

echo [`date +"%Y-%m-%d %H:%M:%S"`]"---------------Running stringtie 1-------------"

cd ~/$job/data/workflow_SE/results/samtools/

for file in $(ls)
	
	do
		stringtie $file -p 50 -v -G $gtf -o ~/$job/data/workflow_SE/results/stringtie/gtfs/string_$file.gtf -A ~/$job/data/workflow_SE/results/stringtie/abundance/stringtie_out_$file.txt
	
	done

echo [`date +"%Y-%m-%d %H:%M:%S"`]"--------Done with stringtie 1 ; Created various GTFs---------"


cd ~/$job/data/workflow_SE/results/stringtie/gtfs/

echo [`date +"%Y-%m-%d %H:%M:%S"`]"--------Running stringtie merge---------"

ls *.gtf > gtf_list

stringtie --merge -G $gtf -o ~/$job/data/workflow_SE/results/stringtie/merge/gtfs/stringstie_merged.gtf gtf_list

echo [`date +"%Y-%m-%d %H:%M:%S"`]"--------Done with stringtie merge---------"

#: << 'END' <COMMENT OUT WHEN USING FEATURE COUNTS>

echo [`date +"%Y-%m-%d %H:%M:%S"`]"--------Running stringtie 2---------"

cd ~/$job/data/workflow_SE/results/samtools/

for file in $(ls)
	
	do
		stringtie $file -p 50 -v -e -G ~/$job/data/workflow_SE/results/stringtie/merge/gtfs/stringstie_merged.gtf -o ~/$job/data/workflow_SE/results/stringtie/gtfs/merge_output_stringtie_$file.gtf -A ~/$job/data/workflow_SE/results/stringtie/merge/abundance/stringtie_merge_out_$file.txt
	
	done

#END <COMMENT OUT WHEN USING FEATURE COUNTS>

echo [`date +"%Y-%m-%d %H:%M:%S"`]"--------Done with all instances of stringtie---------"

#END

##-----------------------------------------##

#: << 'END'

###FeatureCount tool (USING GTF CREATED BY STRINGTIE)

gtf=~/$job/data/workflow_SE/results/stringtie/merge/gtfs/stringstie_merged.gtf #chiranjeevdas gtf=gtf generated by stringtie

echo [`date +"%Y-%m-%d %H:%M:%S"`]"---------------Generating feature counts-------------"

cd ~/$job/data/workflow_SE/results/samtools/

featureCounts -T 50 -t exon -g transcript_id -a $gtf -o counts.txt -M *.bam	#resource limit #chiranjeevdas
	#-T number of threads

echo "---------Done Generating count data-----------"

rsync -a -h -v -r -P -t --remove-source-files *.txt ~/$job/data/workflow_SE/results/featureCounts/
rsync -a -h -v -r -P -t --remove-source-files *.summary ~/$job/data/workflow_SE/results/featureCounts/

echo "--------Results are in featureCount folder of results---------"

echo [`date +"%Y-%m-%d %H:%M:%S"`]"--------Finished with featurecounts-----"

#END

##-----------------------------------------##

#: << 'END'

### DeSeq2 in R

echo [`date +"%Y-%m-%d %H:%M:%S"`]"---------------Calling R for DeSeq2-------------"

cp $deseq2 ~/$job/scripts/deseq2.R
cd ~/$job/scripts/

r deseq2.R

mv res.csv bak_res.csv
echo -n "", > res.csv; cat bak_res.csv >> res.csv #fixes the left shift of column names
rm bak_res.csv 

mv res_PAdj_cutoff.csv bak_res_PAdj_cutoff.csv
echo -n "", > res_PAdj_cutoff.csv; cat bak_res_PAdj_cutoff.csv >> res_PAdj_cutoff.csv #fixes the left shift of column names
rm bak_res_PAdj_cutoff.csv 

echo [`date +"%Y-%m-%d %H:%M:%S"`]"---------------DeSeq2 Complete-------------"

#END

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
