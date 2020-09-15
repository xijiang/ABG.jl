using Statistics

"""
    bgcf()
---
Create scenarios for Theo's `bgcf` program.
"""
function bgcf()
    subtitle("A setup for Theo's bgcf")
    item("Parameters")
    isfile(joinpath(abg_rst, "toy.G")) || simu_qp()
    nlc, nid = countlines.(joinpath.(abg_dat, ["toy.bim", "toy.fam"]))
    open(joinpath(abg_rst, "bgcf.ph"), "w") do io
        for ph in eachline(joinpath(abg_rst, "toy.ph"))
            println(io, "$ph 1 1")
            done
        end
    end
    seed = rand(1:2^31)
    write(joinpath(abg_rst, "dseed.rec"), "$seed\n")
    giv = reshape(reinterpret(Float64, read(joinpath(abg_rst, "toy.giv"))), (nid, nid))
    open(joinpath(abg_rst, "g.giv"), "w") do io
        for i in 1:nid
            for j in 1:nid
                ele = giv[i, j]
                println(io, "$i $j $ele")
            end
        end
    end
    cp(joinpath(abg_dat, "toy.bim"), joinpath(abg_rst, "toy.bim"), force=true)
    cp(joinpath(abg_dat, "toy.bed"), joinpath(abg_rst, "toy.bed"), force=true)
    nchr = 1 
    _π = 0.1
    ve = 0.2
    vsnp = 0.8 * 10 / nlc
    vu = 0.8
    cycle = 1000
    burnin = 100
    nworker = 1
    ithin = 0
    ispeed = 0
    
    inp_str=["&inputvalues",
               "Nanim=$nid",
               "Nchr=$nchr",
               "Nsnptotal=$nlc",
               "bimfil=toy.bim",
               "bedfil=toy.bed",
               "phenfil=bgcf.ph", # line: y weight class; weight=1 → equal weight
               "PI_value=$_π", # for BayesC analysis. π ∈ [0, 1].  π < 0 → estimate π
               "Ve=$ve",
               "Vsnp=$vsnp",  # Vsnp matrix; SNP ∼ N(0, Vsnp) 1% of Vg
               "Vu=$vu", # G=polygenic covariance matrix, 70% of Vg.  All 0 → no polygenic variance. 
               "Ncycle=$cycle",   # Gibbs cycles
               "Nburnin=$burnin", # burnin cycles
               "Nwork=$nworker",  # number of Gibbs sampling threads
               "Ithin=$ithin", # thinning output (prt each ithin sample).  ithin = 0 → no sample output, only smry stats
               "Ispeedup=$ispeed", # reduce computations by factor ispeedup.  ispeedup = 0 → no reduction.
               "/"]
    write(joinpath(abg_rst, "bgcf.inp"), join(inp_str, '\n'), '\n')
    done()
end
