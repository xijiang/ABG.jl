"""
Xijiang's breeding algorithms implemented in Julia and C++.
"""
module ABG
#using SparseArrays, LinearAlgebra, Serialization, Statistics

workdir = pwd()
abg_bin, abg_cpp, abg_dat = begin
    t = splitdir(pathof(ABG))[1]
    l = findlast('/', t) - 1
    d = t[1:l]
    joinpath.(d, ["bin", "cpp", "dat"])
end

global abg_rst = joinpath(pwd(), "sandbox")
set_rst(dir) = (global abg_rst = dir)

beagle = joinpath(abg_bin, "beagle.jar")
plink = joinpath(abg_bin, "plink")

# v0.1
################################################################################
include("styled-messages.jl")   # checked 2020-09-13-11:39
include("final-reports-2-ped.jl")

# v0.2
include("makefile.jl")
include("a-matrix.jl")
include("g-matrix.jl")
include("g-med-stor-m.jl")
include("sort-pedigree.jl")
include("simulation.jl")

# v0.2.1
include("bgcf.jl")              # for Theo's `bgcf` program
include("DMU.jl")               # prepare for the running of DMU

ABG.title("ABG Julia Package")
ABG.message(lpad("Developed by Xijiang Yu @NMBU", 60),
            lpad("MIT license", 60),
            lpad("Copyright © 2020 - " * string(Dates.year(Dates.now())), 60)
            )
ABG.Update()
ABG.message("Current sandbox path is $abg_rst",
            "You can change it by calling 'ABG.set_rst(your_desired_path)'"
            )
end # module
