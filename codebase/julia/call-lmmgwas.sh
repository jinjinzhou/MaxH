#!/bin/sh

#  call-lmmgwas.sh
#  
#
#  Created by Jin Zhou on 9/13/16.
#
Phenotype=fev_std_res
cp copd-example.jl ./split/copd-example_${Phenotype}.jl
sed -ie "s/Phenotype/$Phenotype/g" ./split/copd-example_${Phenotype}.jl

for i in $(seq 22)
do

cp ./split/copd-example_${Phenotype}.jl ./split/copd-example_${Phenotype}_chr${i}.jl
sed -ie "s/chrindex/$i/g" ./split/copd-example_${Phenotype}_chr${i}.jl

#cat > Chr$i.tmp  << EOF
#!/bin/csh

#PBSS -N job${i}
###PBSS -m bea
#PBSS -W group_list=jzhou
#PBSS -q standard
#PBSS -l jobtype=serial
#PBSS -l select=1:ncpus=6:mem=11gb
#PBSS -l pvmem=23gb
#PBSS -l place=pack:shared
#PBSS -l walltime=10:00:00
#PBSS -l cput=10:00:00

### set directory for job execution, ~netid = home directory path
#cd /rsgrps/jzhou/MaxH/MaxH/codebase

### load julia ###
#module load unsupported
#module load markb/julia/0.4.1

### run your executable program with begin and end date and time output
#date
#${gctafolder}gcta64 --bfile ../datasets/CG10kNhwHg19Clean_v2_Mar2013  --chr $i  --make-grm  --out ../datasets/dataCG10kNhwHg19Clean_v2_Mar2013_chr${i}
#date

#EOF
#qsub Chr$i.tmp
echo "./split/copd-example_${Phenotype}_chr${i}.jl"
Julia ./split/copd-example_${Phenotype}_chr${i}.jl
done 

