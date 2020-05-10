module ABG

################################################################################
# export ped_n_map_to_bed, bed_snp_subset, miss_allele_stats, allele_maf_stats
# export hwe_stats, merge_beds, plink_2_vcf, vcf_2_plink, bed_2_map_n_ped
# export plink_filter_snp, plink_filter_id
# export empty_dir
# export title, message, warning, item, done
# export fr2ped
plink = "bin/plink"             # Update this in you own environment

include("styled-messages.jl")
include("final-reports-2-ped.jl")
include("work-flow.jl")
include("misc.jl")
include("plink-cmds.jl")

end # module
