using Statistics
"""
    simu_qp(h2 = 0.8)
---
# Description
Simulation of QTL and phenotypes with the toy dataset comes with this package.
"""
function simu_qp(h2 = 0.8)
    subtitle("Simulation QTL and phenotypes with toy data")
    item("Prepare genotypes")
    isdir(abg_rst) || mkdir(abg_rst)
    nlc, nid = countlines.(joinpath.(abg_dat, ["toy.bim", "toy.fam"]))

    run(`$plink --bfile $abg_dat/toy --recode A --out $abg_rst/toy`)
    run(pipeline("$abg_rst/toy.raw", `$abg_bin/raw2gt $abg_rst/toy.frq`, "$abg_rst/toy.gt"))
    _, twop, Z = read_gt_n_frq("$abg_rst/toy.gt", "$abg_rst/toy.frq")
    done()

    item("QTL, BV, and phenotype simulaiton")
    nq = nlc ÷ 20
    lq = sort(rand(1:nlc, nq))
    eq = randn(nq) .* sqrt(1. / nq)
    bv = Z[lq, :]'eq
    vg = var(bv)
    ve = vg / h2 * (1 - h2)
    er = randn(nid) .* sqrt(ve)
    ph = bv + er
    write(joinpath(abg_rst, "toy.ph"), join(ph, '\n'), '\n')
    write(joinpath(abg_rst, "toy.bv"), join(bv, '\n'), '\n')
    open(joinpath(abg_rst, "toy.qtl"), "w") do io
        for i in 1:nq
            println(io, lq[i], ' ', eq[i])
        end
    end
    done()

    item("G and its inverse")
    G = vanRaden(Z, twop, 1e-4)
    giv = inv(G)
    write(joinpath(abg_rst, "toy.G"), G)
    write(joinpath(abg_rst, "toy.giv"), giv)
    done()
end
