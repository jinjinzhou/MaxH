#!/bin/sh

#  EstimateGVar.sh
#  
#
#  Created by Jin Zhou on 12/13/15.
#
#  This script is used to estimate genetic variances of phenotypes
#  and genetic co-variances between pairs of phenotypes.
#  "fev_std_res"
#  "FVCpp_std_res"
#  "TLCpp_std_res"
#  "fev_fvc_std_res"
#  "logpctEmph_std_res"
#  "logpctEmph_UT_std_res"
#  "logpctEmph_LT_std_res"
#  "logpctGasTrap_std_res"
#  "logUT_LT_std_res"
#  "logPi10_SRWA_std_res"
#   "WallAreaPct_seg_std_res"
#
#  Phenotype file: ../datasets/CG10kNHWPheno4Subtyping.txt
#  Standardized residual file: ../datasets/CG10kNHWRes4Subtyping.txt which adjusted for gender+age+height+packyears+pc1+pc2+pc3
#  Standardized residual file with complete information: ../datasets/CG10kNHWRes4SubtypingNoMissing.txt
#  GRM file: ../datasets/CG10kNhwHg19Clean_v2_Mar2013, note this GRM was recalculated using all genotypes information from COPDGene.
#  Quantitative covariates file: ../datasets/qcovarNHW.txt
#  Categorical covariates file: ../datasets/covarNHW.txt
#
#  Above files were generated using R code: process_phenotype.R
#

#
# If not pre-process covariate files, use the following code to extract qcovar columns
#
#awk -f ../codebase/extract.awk -v cols=FID,IID,Age_Enroll,gender,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 ../datasets/CG10kNHWPhenoPCAs_4plink.txt > tmp.txt
#awk -f ../codebase/extract.awk -v cols=FID,IID,Age_Enroll,gender,PC1,PC2,PC3 ../datasets/CG10kNHWPhenoPCAs_4plink.txt > tmp.txt
#awk -F" " 'BEGIN{OFS=FS;}{t=$i;$i=$j;$j=t;}1' i=2 j=14  tmp.txt > ../datasets/qcovarNHWawk.txt
#rm tmp.txt

#
# call GCTA to estimate genetic variances
#
gctafolder=/rsgrps/jzhou/bin/gcta_1.02/
phenotypefile=CG10kNHWRes4Subtyping.txt
grmfile=CG10kNhwHg19Clean_v2_Mar2013

# create symbolic link to ../datasets/ folder
# ln /rsgrps/jzhou/COPDGene/Phase1_Phenotype/$phenotypefile ../datasets/$phenotypefile
# ln /rsgrps/jzhou/COPDGene/GRM/${grmfile}.grm.gz ../datasets/${grmfile}.grm.gz
# ln /rsgrps/jzhou/COPDGene/GRM/${grmfile}.grm.id ../datasets/${grmfile}.grm.id
#qcovfile=qcovarNHW.txt
#covfile=covarNHW.txt

names=$(head -n 1 ../datasets/$phenotypefile)
namesarray=($names)
len="${#namesarray[@]}"
phenotypes=$((len-2))

for i in $(seq $phenotypes)
do
pi1=$((i+1))
outpi1=${namesarray[$pi1]}
${gctafolder}gcta64 --reml  --grm ../datasets/$grmfile --pheno ../datasets/$phenotypefile --mpheno $i --qcovar ../datasets/$qcovfile --covar ../datasets/$covfile --out $j



ii=$((i+1))
    for j in $(seq $ii $phenotypes)
    do
    pi2=$((j+1))
    outpi2=${namesarray[$pi2]}

cat > job$i"_"$j.tmp  << EOF
#!/bin/csh

#PBS -N job$i"_"$j
###PBS -m bea
#PBS -W group_list=jzhou
#PBS -q standard
#PBS -l jobtype=serial
#PBS -l select=1:ncpus=6:mem=11gb
#PBS -l pvmem=23gb
#PBS -l place=pack:shared
#PBS -l walltime=10:00:00
#PBS -l cput=10:00:00

### set directory for job execution, ~netid = home directory path
cd /rsgrps/jzhou/MaxH/MaxH/codebase

### run your executable program with begin and end date and time output
date
${gctafolder}gcta64 --reml-bivar $i $j --grm ../datasets/$grmfile --pheno ../datasets/$phenotypefile  --out $outpi1$outpi2
date

EOF
	qsub job$i"_"$j.tmp
    done
done



