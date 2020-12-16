#!/bin/bash

##############################################################################################################
################## Storing the path of common address ########################################################
##############################################################################################################


list=/home/cluster/kshattry/Lists

mscan=/home/cluster/kshattry/Matrix-Scan/Arabidopsis_thaliana_PWMs

sequence=/home/cluster/kshattry/Gene_sequences

pseudog=/home/cluster/kshattry/Pseudogenomes/Arabidopsis_thaliana_pseudogenomes

script=/home/cluster/kshattry/Running_scripts

result=/home/cluster/kshattry/Scan_results/Arabidopsis_thaliana

echo "Paths stored"



################################################################################################################
#################### Editing the genelist for any windows-based error/duplications #############################
################################################################################################################

#Converting windows text file (if present) to linux file
tr -d '\15\32' < $list/cumulative > $list/genelist_w1

#Converting lower cases of gene IDs to upper cases
tr '[:lower:]' '[:upper:]' < $list/genelist_w1 > $list/genelist_w2

sort -u $list/genelist_w2 > genelist

rm $list/genelist_w1

rm $list/genelist_w2 

cp $list/genelist $mscan/genelist

###############################################################################################################
################### Retrieving the start and end coordinates of genes #########################################
###############################################################################################################

retrieve-seq  -org Arabidopsis_thaliana.TAIR10.42  -feattype gene -type upstream -format fasta -label id,name -from -1000 -to -1 -i $list/genelist -o $sequence/ref-seq

echo "Upstream regions of the target genes extracted"

cd $pseudog

grep -i ">" $sequence/ref-seq | sed "s/[|>:]/ /g" | awk '{print $1,$17,$18}' | awk '{if($2<$3) print}' > Gene_Coord

echo "Gene Coordinates of the upstream regions extracted"

echo "Files with upstream coordinates transferred to Pseudogenomes"

while read -r line
do

	grep -iw $line $mscan/TF_Information_all_motifs_plus | awk '{if($9=="D") print $4,$6,$7}' >> $mscan/tflist

done < $list/genelist

echo "List of transcription factors extracted"

################################################################################################################
################### Indexing all the pseudogenomes #############################################################
################################################################################################################

#gunzip -v $pseudog/*.gz 

cd $pseudog

#Making a list of all Pseudogenomes
ls *.fasta > eco_list

#Indexing of pseudogenomes using samtools
parallel "samtools faidx {}" < eco_list

################################################################################################################
############################ Extracting the upstream sequences of all genes from all ecotypes ##################
################################################################################################################

for ((i=1;i<=5;i++))
do
	
	awk -v j="AT$i" '{if(substr($1,1,3) == j) print}' $pseudog/Gene_Coord > Gene_Coord_Subset
	while read -r line
	do
		
		a=($line)
        	
		echo $a
        	
		IFS=" " gene=${a[0]}
	  	
		g="${gene}_sequences"
        	
		echo $g
        	
		#Taking the starting coordinate of gene1
        	
		start_c=${a[1]}
        	
		#Taking the ending coordinate of gene1
        	
		end_c=${a[2]}
        	
		#Working on each type of ecotypes
        	
		while read -r line
        	
		do

               	echo "Extracting the upstream sequences of $pseudog/$gene from $line"
                	b=(${line//./' '})
                	acc1=${b[0]}
                	acc=(${acc1//o/' '})
                	name="MPI-GMI|Ath-1001-Genomes|pseudo-genome|${acc[1]}|Chr$i|v3.1"
	         	echo -e "$name\t$start_c\t$end_c" > $pseudog/Coord
			bedtools getfasta -fi $pseudog/$line -bed $pseudog/Coord -fo $pseudog/sequence
			cat $pseudog/sequence >> $sequence/"$g"
	
		done < $pseudog/eco_list

        	sleep 01
       		rm $pseudog/sequence

	done < $pseudog/Gene_Coord_Subset

done



cd $sequence

ls | grep -i "AT" > list_of_sequences

echo "List of sequence made"

###########################################################################################################################
####################### Compilation of all motifs in a file in TRANSFAC format ############################################
###########################################################################################################################

cd $mscan

while read -r line
do
	a=($line)
	a+=".txt"
	#echo $a
	l=($(wc $a))
	if [[ $l == 1 ]]
	then

		sed -i "s/$line//g" tflist

	fi

done < tflist

sed -i '/^[[:space:]]*$/d' $mscan/tflist

echo "Valid list of TFs made"

while read -r line
do

	a=($line)
	motif="${a[0]}.txt"
	echo "AC ${a[0]}" >> $mscan/collected
	echo "XX" >> $mscan/collected
	echo "ID ${a[2]}" >> $mscan/collected
	echo "XX" >> $mscan/collected
	cat $mscan/$motif >> $mscan/collected
	echo "XX" >> $mscan/collected

done < $mscan/tflist

echo "//" >> $mscan/collected

sed -i "s/Pos/PO/g" $mscan/collected

echo "Concatenation of motifs done"

export PATH=$HOME/meme/bin:$HOME/meme/libexec/meme-VERSION:$PATH

cd $sequence

##################################################################################################################################
################### Scanning of sequences by single compiled motif file using FIMO ###############################################
##################################################################################################################################
	
parallel "fasta-get-markov -m 1 -dna {} ./background/{}-bg" < list_of_sequences
parallel "transfac2meme -bg ./background/{}-bg -use_acc $mscan/collected > $mscan/meme_motifs/{}-meme" < list_of_sequences
parallel "fimo -o $result/{}-scan -bfile ./background/{}-bg $mscan/meme_motifs/{}-meme ./{}" < list_of_sequences

rm $mscan/collected

###################################################################################################################################
###################################### Removing duplicated rows from result #######################################################
###################################################################################################################################
cd $result

ls | grep -i "AT" > list_of_scan

while read -r line
do
	
	#cd $line	
	echo $line
	cat $line/fimo.tsv | awk '{print $1,$2,$3,$4,$5,$6,$7,$9}' OFS="\t" | sort -u | grep "M0" > $result/$line/actual_values.tsv
	#cd $result
	
done < list_of_scan

END

##################################################################################################################################
#################################### Replacing the TF ids with name ##############################################################
##################################################################################################################################

cd $result 

ls | grep -i "AT" > list_of_scan

while read -r line
do

	gene=($line)
	echo $gene
	cat $gene/actual_values.tsv | awk '{print $1}' | sort -u > unique
	while read -r line
	do

		tf_id=($line)
		tf=($(grep $tf_id $mscan/tflist))
		sed -i "s/$tf_id/${tf[2]}/g" $gene/actual_values.tsv

	done < unique

done < list_of_scan

rm unique

###################################################################################################################################
##################################### 





