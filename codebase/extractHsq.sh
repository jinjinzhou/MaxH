#!/bin/sh

#  extractHsq.sh
#  
#
#  Created by Jin Zhou on 1/14/16.
#

# define phenotypes and extract phenotype names
# first two colums are FID and IID
phenotypefile=CG10kNHWRes4Subtyping.txt
names=$(head -n 1 ../datasets/$phenotypefile)
namesarray=($names)
len="${#namesarray[@]}"
phenotypes=$((len-2))


cat > estHsq.txt << EOF
EOF

for i in $(seq $phenotypes)
do
echo $i
# get the phenotypes name
pi1=$((i+1))
outpi1=${namesarray[$pi1]}
echo $outpi1
hsqfile=${outpi1}.hsq
gvar=$(grep -r 'V(G)'  $hsqfile)
n=$(grep -A 1 'Pval' $hsqfile | grep  'n')
echo ${outpi1}" "$gvar" "$n >> estHsq.txt
ii=$((i+1))
    for j in $(seq $ii $phenotypes)
    do
    echo $j
    pi2=$((j+1))
    outpi2=${namesarray[$pi2]}
    hsqfile=$outpi1$outpi2".hsq"
    gvar=$(grep -r 'C(G)_tr12'  $hsqfile)
    n=$(grep -A 1 'Pval' $hsqfile | grep  'n')
    echo $outpi1$outpi2" "$gvar" "$n >> estHsq.txt
    done
done


