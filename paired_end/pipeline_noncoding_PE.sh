#!/bin/bash

#SBATCH --job-name=tomato
#SBATCH --nodes=1
#SBATCH --ntasks=56

#################################################################
#################     USER INSTRUCTIONS     #####################
#################################################################

# 1. Modify the "job-name" above (OPTIONAL)
# 2. In the "DECLARATIONS/PATHS" section below, put a project name.
# 3. In the "DECLARATIONS/PATHS" section below, put a unique job name. NB. This is specific to each run of the pipeline and must be edited every time you run it.
# 4. Populate the 'sra', 'fatstq_f', 'ref_g', 'gtf_f' & 'deseq2' paths. Ensure the paths/files exist.
# 5. Save the edited script as a copy and run with sbatch <script name>

# IMPORTANT: The control files of the data set must have the work "control" in the filename
# This is used for DESeq2. e.g. control101.fastq.gz, 1_control.sra, 1controlb.fastq.gz

#################################################################
#################     USER SCRIPT START     #####################
#################################################################

start_time=`date +%s`
echo "start time:" $start_time

#################################################################
#################   DECLARATIONS / PATHS   ######################
#################################################################
proj="test" #manually set job name here before running, one project can have multiole jobs
job="tomato" #manually set job name here before running

#sra="/home/cluster/akash/projects/rsolani_nopaper/data/srarawdata" #sra files (source) directory (no trailing /)
fastq_f="/home/cluster/drishtee/data/crops/TOMATO/cold_drought/ena_files/cold_drought/" #fastq files (source) directory (no trailing /)
ref_g="/home/cluster/drishtee/data/reference_genomes/tomato_ref_gen.fa" #full path of your (source) reference genome file ending in .fa OR .fasta
gtf_f="/home/cluster/drishtee/data/gtf_files/Solanum_lycopersicum.SL3.0.49.gtf" #full path of (source) gtf file ending in .gtf
#deseq2="/path/to/deseq2/R/script.R" #full path of the DeSeq2 R script (.R)

#################################################################
################ Start of Pipline workflow_PE ###################
#################################################################
echo [`date +"%Y-%m-%d %H:%M:%S"`] "----STARTING THE PIPELINE------"

#: <<'END'
echo [`date +"%Y-%m-%d %H:%M:%S"`] "----------Creating Directory hierachy------"

mkdir -p ~/projects/$proj/pipeline_result/$job/scripts/
mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastqc/
mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastq_files/
mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastp/
mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastp_reports/
mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/featureCounts/
mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/hisat2/
mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/samtools/
mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/stringtie/gtfs/
mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/stringtie/abundance/
mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/stringtie/merge/gtfs/
mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/stringtie/merge/abundance/
mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/stringtie/ctabs/
mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/stringtie/all_ctabs/
mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/degs/stringtie
mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/degs/featurecounts
mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/deseq/stringtie
mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/deseq/featurecounts
mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/reference_genome/
mkdir -p ~/projects/$proj/pipeline_result/$job/data/workflow_PE/gtf/

: << 'END'
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
END

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

echo [`date +"%Y-%m-%d %H:%M:%S"`] "---------------Generating feature counts with string tie-------------"

#: << 'END'

cd ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/samtools/

for file in $(ls *.bam)
	do
		stringtie $file -p 50 -v -G $gtf_f -o ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/stringtie/gtfs/string_$file.gtf -A ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/stringtie/abundance/stringtie_out_$file.txt
	done

echo [`date +"%Y-%m-%d %H:%M:%S"`] "---------------merging stringtie gtfs-------------"

cd ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/stringtie/gtfs/

ls *.gtf > gtf_list

stringtie --merge -G $gtf_f -o ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/stringtie/merge/gtfs/stringstie_merged.gtf gtf_list

echo [`date +"%Y-%m-%d %H:%M:%S"`] "---------------Generating ctabs -------------"

cd ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/samtools/

for file in $(ls)

	do
		stringtie -p 50 -v -e -b ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/stringtie/ctabs/$file/ -G ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/stringtie/merge/gtfs/stringstie_merged.gtf $file -o ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/stringtie/gtfs/merge_output_stringtie_$file.gtf -A ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/stringtie/merge/abundance/stringtie_merge_out_$file.txt 
	done

echo [`date +"%Y-%m-%d %H:%M:%S"`] "---------------Consolidating ctabs-------------"

cd ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/stringtie/

find ./ -name *t_data.ctab > temp
#finds all paths to files with "t_data.ctab" in name and writes them to a "temp" file
echo "found `cat temp | wc -l` ctabs at:"
cat temp
echo "copying to consolidated location"
#loop to read the paths from temp file one by one and store in "x"
for x in $(cat temp)
do
cp -v $x ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/stringtie/all_ctabs/`echo $x | awk -F "/" '{print $(NF-1)}' | sed s/.bam//`.ctab
#copies each file (using path stored in "x") to the folder "all_ctabs" and renames it according to source file name.
done #loop end

rm temp

echo [`date +"%Y-%m-%d %H:%M:%S"`] "---------------stringtie complete-------------"

#END

#: << 'END'

echo [`date +"%Y-%m-%d %H:%M:%S"`] "---------------creating r script to read ctabs and give DEGs-------------"

cat > ~/projects/$proj/pipeline_result/$job/scripts/ctab_to_deseq2.r<< EOF
#!/usr/bin/r
library("stringr")
library("DESeq2")
library("ggplot2")
library("dplyr")
library("readr")
library("tximport")
library("gplots")
setwd("/projects/$proj/pipeline_result/$job/data/workflow_PE/results/stringtie/all_ctabs") #path to ctabs
getwd()
list.files()
#Add all .ctab files in a vector
files <- list.files(getwd(), pattern=".ctab", all.files=FALSE, full.names=FALSE)
#Reorder control and sample files to match dataframe structure
f1 <- grep("control",files)
f1 <- files[f1]
f2 <- files[!(files %in% f1)]
f3 <- union(f2,f1)
files <- f3
# count the sample and control files
no_cont <-length(grep("control",files,ignore.case=TRUE))
no_samp <- (length(files) - no_cont) 
#Read the ctab files, and store in txi
tmp <- read_tsv(files[1])
tx2gene <- tmp[, c("t_name", "gene_id")]
txi <- tximport(files=files, type = "stringtie", tx2gene = tx2gene) 
# Define conditions for the samples
sampleTable <- data.frame(condition = factor(c(rep("Sample",no_samp), rep("Control",no_cont))))
sample_names <- sapply(strsplit(files, "[.]"), `[`, 1)
row.names(sampleTable) <- sample_names
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
dds <- dds[rowSums(counts(dds)) > 10,]
dds <-DESeq(dds)
vst <- vst(dds, blind=FALSE)
# Plot PCA plot
svg("stie_pca_vst.svg")
plotPCA(vst, intgroup="condition", ntop=nrow(counts(dds)))
dev.off()  
# Explore PCA plot
svg("stie_pca_vst_2.svg")
a <- DESeq2::plotPCA(vst, intgroup="condition")
a + geom_label(aes(label = sampleTable$condition),)
nudge <- position_nudge(y = 1)
a + geom_label(aes(label = sampleTable$condition), position = nudge)
a + geom_text(aes(label = sampleTable$condition), position = nudge, size=3 )
dev.off()
svg("stie_boxplot_vst.svg")
boxplot(assay(vst), outline = FALSE, main = "Boxplot based on vst transformation", font.main= 2, font.axis=0.5, font.lab=2, col=27, col.axis=2, cex=0.5, ylim=c(-10,11))
dev.off() 
svg("stie_boxplot_samples.svg")
boxplot(assay(vst), col= c("Red"), pch=".",
vertical=TRUE, cex.axis=0.5, main = "Boxplot of samples using vst method",
las=2, ylab="assay(vst)", xlab="Samples", ylim=c(-10,30),font.main= 5, font.axis=0.5, font.lab=2 )
dev.off() 
# Plot correlation heatmap
cU <-cor( as.matrix(assay(vst)))
cols <- c("dodgerblue3", "firebrick3")[sampleTable$condition]
svg("stie_heatmap.svg")
heatmap.2(cU, symm=TRUE, col= colorRampPalette(c("darkblue","white"))(100),
            labCol=colnames(cU), labRow=colnames(cU),
            distfun=function(c) as.dist(1 - c),
            trace="none",
            Colv=TRUE, cexRow=0.9, cexCol=0.9, key=F,
            font=2,
            RowSideColors=cols, ColSideColors=cols)
dev.off() 
#dispersion plot
svg("stie_dispersion.svg")
plotDispEsts(dds)
dev.off() 
res <- results(dds, contrast=c("condition","Sample","Control" ))
summary(res)
grp.mean <- sapply(levels(dds$condition),
                          function(lvl)
                            rowMeans(counts(dds,normalized=TRUE)[,dds$condition== lvl]))
norm.counts <- counts(dds, normalized=TRUE)
all <- data.frame(res, assay(vst))
nrow(all)
write.table(all, file="stie_degs.csv",sep=",")
padj_cutoff <- all[all$padj <= 0.05,]
write.table(padj_cutoff, file="stie_degs_padj.csv",sep=",")
write.table(assay(vst), file="stie_vst_table.csv",sep=",")
EOF

echo [`date +"%Y-%m-%d %H:%M:%S"`] "---------------executing the ctab to deseq2 script-------------"

cd ~/projects/$proj/pipeline_result/$job/scripts/

r ctab_to_deseq2.r

sleep 15

mv -v ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/stringtie/all_ctabs/stie_* ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/deseq/stringtie/

#END

echo [`date +"%Y-%m-%d %H:%M:%S"`] "---------------Fixing column leftshift---{{{{{{{{{{{{{{{{{{{{{{{TO BE DONE}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}--"

cd ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/deseq/stringtie/
mv stie_degs.csv bak_res.csv
echo -n "", > stie_degs.csv; cat bak_res.csv >> stie_degs.csv #fixes the left shift of column names
rm bak_res.csv 

mv stie_degs_padj.csv bak_res_PAdj_cutoff.csv
echo -n "", > stie_degs_padj.csv; cat bak_res_PAdj_cutoff.csv >> stie_degs_padj.csv #fixes the left shift of column names
rm bak_res_PAdj_cutoff.csv 


##-----------------------------------------------##

echo [`date +"%Y-%m-%d %H:%M:%S"`]"---------------Cleaning Up-------------"

echo [`date +"%Y-%m-%d %H:%M:%S"`]"---------------Deleting SAM files-------------"
rm -v ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/hisat2/*
echo [`date +"%Y-%m-%d %H:%M:%S"`]"---------------Deleting Trimmed FastQ Files-------------"
rm -v ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastp/*
echo [`date +"%Y-%m-%d %H:%M:%S"`]"---------------Deleting FastQ files-------------"
rm -v ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/fastq_files/*
echo [`date +"%Y-%m-%d %H:%M:%S"`]"---------------Deleting redundant bam files-------------"
rm -v ~/projects/$proj/pipeline_result/$job/data/workflow_PE/results/stringtie/ctabs/*.bam
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
