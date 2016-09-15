using SnpArrays
using DataFrames

cd("$(homedir())/Documents/Research/MaxH/MaxH/codebase/")
@time copd = SnpArray("../datasets/CG10kNhwHg19Clean_v2_Mar2013")
srand(123)
@time Φgrm = grm(copd)
using JLD
save("./mm/cg10k.grm.jld","Φgrm",Φgrm)
