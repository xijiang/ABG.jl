"""
Xijiang's breeding algorithms implemented in Julia and C++.
"""
module ABG
using SparseArrays, LinearAlgebra, Serialization, Statistics

abgDir = begin
    t = splitdir(pathof(ABG))[1]
    l = findlast('/', t)
    SubString(t, 1:l)
end
abgBin = joinpath(abgDir, "bin")
abgCpp = joinpath(abgDir, "cpp")

# v0.2
include("makefile.jl")
include("a-matrix.jl")
include("g-matrix.jl")
include("sort-pedigree.jl")

# v0.1
################################################################################
# export ped_n_map_to_bed, bed_snp_subset, miss_allele_stats, allele_maf_stats
# export hwe_stats, merge_beds, plink_2_vcf, vcf_2_plink, bed_2_map_n_ped
# export plink_filter_snp, plink_filter_id
# export empty_dir
# export title, message, warning, item, done
# export fr2ped
include("styled-messages.jl")
include("final-reports-2-ped.jl")
include("work-flow.jl")
include("misc.jl")
include("plink-cmds.jl")

end # module
