# dm6 pair-end read depth calculate
# !/bin/bash

end=$1
length=$2
genome_len=$3

function printHelp(){
    echo "Version: 1.0"
    echo
    echo
    echo "Required Arguments:"
    echo
    echo " $ 1                     paired-end is 2 and single is 1"
    echo " $ 2                     reads length is default 150 bp or others"
    echo " $ 3                     dm6 genome length is 146638899"    
    echo

}


echo "calculating the sequencing depth ~~~"
echo Start time: `date`
echo "practice makes perfect!"
echo "contact with the author: jl_chen@csu.edu.cn"

ls *gz > sample.txt

for i in $(cat sample.txt);do

    zcat $i |wc |awk '{print ("'${end}'"*"'${length}'"*($1/4))}' > seq_len.txt

done

for i in $(cat seq_len.txt );do
    
    cat $i | while read line ;do awk '{prnit $1/"'${genome_len}'" }' > read_depth.txt  #双端测序要*2，读长150bp要*150,果蝇基因组长度约为146638899 

done

echo End time: `date`
