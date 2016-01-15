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
phenotypefile=CG10kNHWRes4Subtyping.txt
grmfile=copd_nhw_filt
qcovfile=qcovarNHW.txt
covfile=covarNHW.txt

names=$(head -n 1 ../datasets/$phenotypefile)
namesarray=($names)
len="${#namesarray[@]}"
phenotypes=$((len-2))

for i in $(seq $phenotypes)
do
pi1=$((i+2))
outpi1=${namesarray[$pi1]}
#/Users/jzhou/Documents/Bin/gcta_1.02/gcta_mac --reml  --grm ../datasets/$grmfile --pheno ../datasets/$phenotypefile --mpheno $i --qcovar ../datasets/$qcovfile --covar ../datasets/$covfile --out $j
#/Users/jzhou/Documents/Bin/gcta_1.02/gcta_mac --reml  --grm ../datasets/$grmfile --pheno ../datasets/$phenotypefile --mpheno $i --out $outpi1
ii=$((i+1))
    for j in $(seq $ii $phenotypes)
    do
    pi2=$((j+2))
    outpi2=${namesarray[$pi2]}
    /Users/jzhou/Documents/Bin/gcta_1.02/gcta_mac --reml-bivar $i $j --grm ../datasets/$grmfile --pheno ../datasets/$phenotypefile  --out $outpi1$outpi2
    done
done



