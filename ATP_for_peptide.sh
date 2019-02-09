#!/bin/bash

clear
DIR=`pwd`
echo "Program written by Abhijeet Singh (abhijeetsingh.aau@gmail.com)
"
if [ "$#" -gt 0 ]; then
	filename=$1
    file=$1
else
echo "Enter your multifasta file name"
read -p 'Filename: ' file

fi
    if [ -z "$file" ]
        then
            exit
    fi
    if [ ! -f $file ]; then
    echo "
ERROR: File \"$file\" not found in \"$DIR\"  
"
    exit 0
    fi
echo "
Processing: $file
"
mkdir -p analysis 
cp ${file} ./analysis/
cd analysis
sed -e 's/ /_/g;s/-/_/g;s/,/_/g' ${file} > ${file}.tmp
cat ${file}.tmp |\
awk '/^>/ {if(N>0) printf("\n"); printf("%s\n",$0);++N;next;} { printf("%s",$0);} END {printf("\n");}'|\
split -l 2 --additional-suffix=.fa - Peptide_
ls Peptide_*.fa > fasta_list
for i in $(cat fasta_list)
    do 
        awk 'BEGIN{RS=">"}NR>1{gsub("\n","\t");gsub("\n","");print RS$0}' ${i} > ${i}_tab
        awk '{ print $2 }' ${i}_tab > ${i}_col2
        cat ${i}_col2 | wc -c > ${i}_len
        var1=`cat ${i}_len` 
        var2=2
        echo  $(($var1 - $var2)) >> peptide_bonds
        
    done
echo "Calculation ATP equivalant for peptide bonding"
awk '{ print $1 * -4 }' peptide_bonds > ATPeq
sed -e 's/^/and?(ATP equivalant for peptide_bonds)?=?/g' ATPeq > ATPeq_2
#https://mcb.berkeley.edu/labs/garcia/sites/mcb.berkeley.edu.labs.garcia/files/Teaching/2017-MCB137/Stouthamer1973.pdf
echo "A	Alanine	+1
C	Cysteine	-3	
D	Aspartic_Acid	0
E	Glutamic_Acid	+1
F	Phenylalanine	-2
G	Glycine	0
H	Histidine	-7	
I	Isoleucine	-1
K	Lysine	0
L	Leucine	+3	
M	Methionine	-4	
N	Asparagine	-2
P	Proline	0
Q	Glutamine	0
R	Arginine	-3
S	Serine	0
T	Threonine	-2
V	Valine	+2
W	Tryptophan	-5
Y	Tyrosine	-2" > list
awk '{ print $1 }' list > AA
awk '{ print $3 }' list > AA_E
echo "Calculating ATP consumed or formed in the aminoacid biosynthesis"
for l in $(ls Peptide_*_col2)
do
    for A in $(cat AA)
    do
        cat $l | grep -o "$A" | wc -l
        
    done > ${l}_values
done 
echo "Preparing final results"
ls Peptide_*_values > value_list
for v in $(cat value_list)
do
    paste ${v} AA_E > ${v}_AAE
    awk '{ print $1 * $2 }' ${v}_AAE > ${v}_AAE_prod
    awk '{ sum += $1; }  END { print sum; }' ${v}_AAE_prod > ${v}_AAE_prod_prod_sum
    cat ${v}_AAE_prod_prod_sum >> total
done
sed 's/^/*/g' total > total_ATP
grep -e ">" ${file} | cut -d " " -f1 > header
paste header total_ATP > energy
paste energy ATPeq_2 > energy_2
sed -e 's_*_-->?(moles ATP/mole precursor)?=?_g;s_?_\t_g' energy_2 > energy_3
sed -e 's/?/\t/g' energy_3 > energy_4
echo -e "Net ATP consumed(-) or formed(+) in the synthesis of total aminoacid precursors 
and ATP equivalant consumed(-) in the formation of peptide bonds for:
(approximately!)
\n$(cat energy_4)" > ${file}_energy_of_formation.txt
clear
cat ${file}_energy_of_formation.txt
cp ${file}_energy_of_formation.txt ../
cd ..
rm -fr analysis
echo "
Results are approximation!
"
echo "                 --------------------------------"
echo "Result file --->> ${file}_energy_of_formation.txt"
echo "                 --------------------------------"
