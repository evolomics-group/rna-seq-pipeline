#!/bin/bash

# rna seq pipeline prerequisites install script (bash)
# this script has been tested to work with Ubuntu 20.04.2 LTS
# other OS versions may require manual installation of the applications

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
echo [`date +"%Y-%m-%d %H:%M:%S"`] "installing apt packages now"
echo "##############################################################"

sudo apt install pigz
sleep 05
sudo apt install fastp
sleep 05
sudo apt install fastqc
sleep 05
sudo apt install hisat2
sleep 05
sudo apt install samtools
sleep 05
sudo apt install stringtie
sleep 05
sudo apt install subread
sleep 05
sudo apt install r-base
sleep 05
sudo apt install littler
sleep 05
sudo apt install r-cran-biocmanager
sleep 05
sudo apt install r-cran-stringr
sleep 05
sudo apt install r-bioc-deseq2
sleep 05
sudo apt install r-cran-ggplot2
sleep 05
sudo apt install r-cran-dplyr
sleep 05
sudo apt install r-cran-readr
sleep 05
sudo apt install r-bioc-tximport
sleep 05
subo apt install r-cran-gplots
sleep 05

echo "##############################################################"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "apt packages installation complete"
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
