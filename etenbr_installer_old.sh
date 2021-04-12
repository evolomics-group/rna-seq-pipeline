#!/bin/bash

# rna seq pipeline prerequisites install script for older systems (bash)

echo "NOTE: root access is necessary for the script"
echo "NOTE: for cluster environment, the script must be run on each node"
echo "NOTE: access to the world wide web is necessary for all installations to complete"
echo "your netork must have access to apt, cran, bioconductor and python repositories"

sleep 05

echo "IMPORTANT: many tools require user input in the install stage, please be present to select relevant options"

sleep 10

mkdir etenbr_tools
cd etenbr_tools

echo "##############################################################"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "updating apt repositories"
echo "##############################################################"

sudo apt update
sleep 05

echo "##############################################################"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "apt repositories updated"
echo "##############################################################"

echo "##############################################################"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "downloading sra-toolkit now"
echo "##############################################################"

wget --output-document sratoolkit.tar.gz http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
sleep 05

echo "##############################################################"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "installing sra-toolkit now"
echo "##############################################################"

tar -vxzf sratoolkit.tar.gz
sleep 05

sratoolkit=`ls -d sra*/`
export PATH=$PATH:$PWD/$sratoolkit'bin'
echo "export PATH=\$PATH:$PWD/$sratoolkit"'bin' >> ~/.bashrc
vdb-config -i

echo "sratoolkit installed in `which fastq-dump`"

echo "##############################################################"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "sra-toolkit installation complete"
echo "##############################################################"

echo "##############################################################"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "installing pigz now"
echo "##############################################################"

sudo apt install pigz
sleep 05


echo "##############################################################"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "pigz installation complete"
echo "##############################################################"

echo "##############################################################"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "installing fastp now"
echo "##############################################################"

sudo apt install fastp
sleep 05


echo "##############################################################"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "fastp installation complete"
echo "##############################################################"

echo "##############################################################"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "installing fastqc now"
echo "##############################################################"

sudo apt install fastqc
sleep 05

echo "##############################################################"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "fastqc installation complete"
echo "##############################################################"

echo "##############################################################"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "installing hisat2 now"
echo "##############################################################"

sudo apt install hisat2
sleep 05

echo "##############################################################"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "  installation complete"
echo "##############################################################"

echo "##############################################################"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "installing samtools now"
echo "##############################################################"

sudo apt install samtools
sleep 05

echo "##############################################################"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "samtools installation complete"
echo "##############################################################"

echo "##############################################################"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "installing stringtie now"
echo "##############################################################"

sudo apt install stringtie
sleep 05

echo "##############################################################"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "stringtie installation complete"
echo "##############################################################"

echo "##############################################################"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "installing featurecounts (subread) now"
echo "##############################################################"

sudo apt install subread
sleep 05

echo "##############################################################"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "subread installation complete"
echo "##############################################################"

echo "##############################################################"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "installing R now"
echo "##############################################################"

sudo apt install r-base
sleep 05

sudo apt install littler
sleep 05

echo "##############################################################"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "R installation complete"
echo "##############################################################"

echo "##############################################################"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "creating R script for R packages installation"
echo "##############################################################"

cat > etenbr_installer.r<< EOF
#!/usr/bin/r
##	IF YOUR SYSTEM IS BEHIND A PROXY, MODIFY THE FOLLOWING AND UNCOMMENT THEM
#	Sys.setenv(http_proxy="http://user:passwd@url:port")
#	Sys.setenv(http_proxy="https://user:passwd@url:port")

update.packages(ask = FALSE)
install.packages("readr")
install.packages("stringr")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("ashr")
install.packages("gplots")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
library("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("tximport")
BiocManager::install("apeglm")
BiocManager::install("vsn")

EOF

sleep 05

echo "##############################################################"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "R script created and written to etenbr_installer.R"
echo "##############################################################"


echo "##############################################################"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "executing R script for package installation"
echo "##############################################################"

sudo r etenbr_installer.r
sleep 05

echo "##############################################################"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "R packages installation complete"
echo "##############################################################"

echo "##############################################################"
echo "##############################################################"
echo "##############################################################"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "installer script completed"
echo "##############################################################"
echo "##############################################################"
echo "##############################################################"

sleep 05
exit
