using LmmGwas, SnpArrays, VarianceComponentModels
using DataFrames
using JLD

cd("$(homedir())/Documents/Research/MaxH/MaxH/codebase/")
# load pre calcuated grm
@load "./mm/cg10k.grm.jld"

# load datasets
@time copd = SnpArray("../datasets/CG10kNhwHg19Clean_v2_Mar2013")
@time copd_snpdata = SnpData("../datasets/CG10kNhwHg19Clean_v2_Mar2013")
people, snps = size(copd)
@time maf, minor_allele, missings_by_snp, missings_by_person = summarize(copd)

basepairs = copd_snpdata.basepairs
chromosome = map(x->parse(Int64,x),copd_snpdata.chromosome)
snpnames = copd_snpdata.snpid

cg10k_trait = readtable("../datasets/CG10kNHWRes4SubtypingImputed.txt";
    separator = ' ',
    header=true,
    eltypes = [UTF8String; UTF8String; Float64; Float64; Float64; Float64; Float64;
               Float64; Float64; Float64; Float64; Float64; Float64; Float64; Float64]);
phenonames = names(cg10k_trait);
#phenotype will be replaced by outside input.
phenoindex = findfirst(phenonames,:Phenotype);
Y = convert(Vector{Float64}, cg10k_trait[:,phenoindex]);

#chrindex will be replace by outside input.
copd_chrindex = copd[:,chromosome.==chrindex];
pvalue_model_0 = lmmGWAS(Y,ones(people),copd_chrindex,Φgrm);

## output results
filename="./mm/Phenotype_lmmgwas_chrindex.txt";
out = DataFrame(snpnames=snpnames[chromosome.==chrindex],chromosome=chromosome[chromosome.==chrindex],basepairs=basepairs[chromosome.==chrindex],pvalue=pvalue_model_0)
writetable("$filename", out);


### now analyze data using diagolized grm ###
(degree_copd, neighbor_copd) = GRMtoGraph(Φgrm,0.008)
(component_copd, next_copd) = findconnect(degree_copd, neighbor_copd);
(sortindex_copd,sortedGRM_copd) = diagonalize_grm(Φgrm,component_copd,next_copd)
length(unique(component_copd))

pvalue_model_1 = lmmGWAS(Y[sortindex_copd],ones(people),copd_chrindex,sortedGRM_copd)
## output results
filename="./mm/Phenotype_gwaslmm_diagonalized_chrindex.txt";
out = DataFrame(snpnames=snpnames[chromosme.==chridex],chromosome=chromosome[chromosme.==chridex],basepairs=basepairs[chromosme.==chridex],pvalue=pvalue_model_1)
writedlm("$filename", out);
