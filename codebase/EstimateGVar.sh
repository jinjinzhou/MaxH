#!/bin/sh

#  EstimateGVar.sh
#  
#
#  Created by Jin Zhou on 12/13/15.
#
#  This script is used to estimate genetic variances of phenotypes
#  and genetic co-variances between pairs of phenotypes.
#
#  "logpctEmph_UpperLobes",
#   "logpctEmph_LowerLobes",
#   "logpctEmph_UL_LL_ratio",
#   "logpctEmph_Slicer",
#   "logSlicer_15pctIn_Total",
#   "logpctGasTrap_Slicer",
#   "logpctEmph_UpperThird_Slicer",
#   "logpctEmph_LowerThird_Slicer",
#   "logPi10_SRWA",
#   "WallAreaPct_seg",
#   "TLCpp_race_adjusted",
#   "logFRCpp_race_adjusted",
#   "FEV1pp_utah",
#   "FVCpp_utah",
#   "FEV1_FVC_utah",
#   "BDR_pct_FEV1",
#   "BDR_pct_FVC"
#
#  Phenotype file: ../datasets/CG10kNHWPheno4Subtyping.txt
#  Standardized residual file: ../datasets/CG10kNHWRes4Subtyping.txt which adjusted for gender+age+agesquare+age*gender+agesquare*gender+pc1+pc2+pc3 
#  GRM file: ../datasets/copd_nhw_filt, note this GRM was calcuated and filted based on Jin Zhou's AJCRM heritability paper.
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
grmfile=copd_nhw_filt

# create symbolic link to ../datasets/ folder
ln /rsgrps/jzhou/COPDGene/Phase1_Phenotype/$phenotypefile ../datasets/$phenotypefile
ln /rsgrps/jzhou/COPDGene/GRM/${grmfile}.grm.gz ../datasets/${grmfile}.grm.gz
ln /rsgrps/jzhou/COPDGene/GRM/${grmfile}.grm.id ../datasets/${grmfile}.grm.id
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
#/Users/jzhou/Documents/Bin/gcta_1.02/gcta_mac --reml  --grm ../datasets/$grmfile --pheno ../datasets/$phenotypefile --mpheno $i --qcovar ../datasets/$qcovfile --covar ../datasets/$covfile --out $j

cat > job$i.tmp  << EOF
#!/bin/csh

#PBS -N job${i}
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
${gctafolder}gcta64 --reml  --grm ../datasets/$grmfile --pheno ../datasets/$phenotypefile --mpheno $i --out $outpi1
date

EOF
qsub job$i.tmp

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



