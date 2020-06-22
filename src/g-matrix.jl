"""
    function read_gt_n_frq(gt::AbstractString, frq::AbstractString)
---
Read genotypes and their frequencies from file `gt` and `frq`.
Note gt can be prepared like below:
```
plink --species --bfile source --recode A --out target
cat target.raw | pathto/raw2gt target.frq > target.gt
```
See also the function `prepare_gt_from_plink_files`.
"""
function read_gt_n_frq(gt::AbstractString, frq::AbstractString)
    title("Read genotypes and frequencies")
    item("Frequency related")
    # Frequency related
    p = Float64[]
    for x in eachline(frq)
        push!(p, parse(Float64, x))
    end
    twop = p.*2.
    done()

    item("Genotypes")
    nlc = length(readline(gt))
    fsz = stat(gt).size
    lsz = nlc + 1               # assuming Unix text format
    nid = Int(fsz/lsz)
    if(nid*lsz ≠ fsz)
        warning("Not a square file")
        return
    end
    message("NID: $nid;\tN_loci: $nlc")
    z = Array{Float64, 2}(undef, nlc, nid) # one column per ID
    fgt = open(gt, "r")
    for j in 1:nid
        t = read(fgt, lsz)
        for i in 1:nlc
            z[i, j] = Float64(t[i] - 0x30) # 0x30 == '0'
        end
    end
    close(fgt)
    done()
    return p, twop, z
end

"""
    vanRaden(Z::Array{Float64, 2})
---
This function calculate **`G`** with van Raden 2008 on genotype **`Z`**.
method I.  Note **`Z`** here is ID column majored.

`G = Z`Z/2Σp_i(1-p_i)`

where **`Z`** contains the marker genotypes for all animals at all loci, corrected for 
the allele frequency per locus, and `p_i` is the frequency of the allele at locus 
*i* coded with the highest integer value.  **`Z`** is derived from the genotypes of the 
animals by subtracting 2 times the allele frequency, that is `2p_i`, from matrix 
**`X`**, which specifies the marker genotypes for each individual as 0, 1 or 2. Values 
for `p_i` are calculated from the data (default), or can be specified in a file by 
the user.

Pass copy(GT) to this function to avoid GT matrix modification.
"""
function vanRaden(Z::Array{Float64, 2}, twop::Array{Float64, 1})
    title("Calculate G with vanRaden method I")
    Z .-= twop
    s2pq = (1 .- .5 .* twop)'twop
    r2pq = 1. / s2pq
    G = Z'Z .* r2pq
end

"""
    vanRaden2(W::Array{Float64, 2})
---
Calculate **`G`** with van Raden method II.
"""
function vanRadenII(W::Array{Float64, 2})
    title("Calculate G with van Raden method II")
    twop = mean(W, dims=2)
    topq = sqrt.(twop .* (1 .- .5 .* twop))
    W .-= twop
    W ./= topq
    G  = W'W ./ length(twop)
end

"""
    yangG(W::Array{Float64, 2})
---
Calculate **`G`** with Yang Jian's method.
"""
function yangG(W::Array{Float64, 2})
    title("Calculate G with Yang's method")

    message("Will address this later")
    # twop = mean(W, dims=2)
    # topq = sqrt.(twop .* (1 .- .5 .* twop))
    # nid = size(gt)[2]
    # for i in 1:nid
    #     
    # 
    # # Calculate the diagonals first
    # for i in 1:lengt
    # W .-= twop
    # W ./= topq
    # # Above is just like van Raden method II
    # for i in 1:length(twop)
end

"""
    D_matrix(M::Array{Float64, 2})
---
Calculate **`D`** (dominance relationship) matrix with genotypes given.
"""
function D_matrix(M::Array{Float64, 2})
    title("Dominance relationship matrix")
    p = mean(M, dims=2).* .5
    q = 1 .- p
    op = [vec(-2 .* p .* p) vec(2 .*p .* q)  vec(-2 .*q .* q)]
    for j in 1:size(M)[2]
        for i in 1:size(M)[1]
            k = Int(M[i, j]) + 1
            M[i, j] = op[i, k]
        end
    end
    s2pq = sum((p .* q .* 2) .^2)
    D = M'M ./s2pq
end

"""
    function plink_2_vr_I(file::AbstractString)
---
Given the file stem name `file`, this funciton calculate `G` with vanRaden method I.
Note, data are stored in `file.ped`, `file.bim`, `file.fam`.
"""
function plink_2_vr_I(file::AbstractString, species::AbstractString="cow")
    title("Calculate G from plink data with vanRaden method I")
    message("   data: $file.{bed,bim,fam}\nSpecies: $species")
    if Sys.which("plink") == Nothing
        warning("Command plink can't be find in default paths")
        return
    end
    isfile(joinpath(abgBin, "raw2gt")) || make()
    tmp = joinpath(workdir, "tmp")
    isdir(tmp) || mkdir(tmp)
    plink_012(file, joinpath(tmp, "plink"), species)
    _ = run(pipeline(joinpath(tmp, "plink.raw"), `$abgBin/raw2gt $tmp/plk.frq`, "$tmp/plk.gt"))
    _, twop, z = read_gt_n_frq(joinpath(tmp, "plk.gt"), joinpath(tmp, "plk.frq"))
    G = vanRaden(z, vec(twop))
    done()
    return G
end

"""
Note: function names start with a `_` is just for testing.
"""
function _test_grm()
    nlc = 1000
    nid = 10
    M = rand(0:2, nlc, nid)
    p = mean(M, dims=2) .* .5
    twop = 2p
    s2pq = sum((1 .- p) .* twop)
    r2pq = 1. / s2pq
    Z = M .- twop
    if Sys.which("calc_grm") ≠ Nothing
        inp = ["$nlc", "genotypes.txt", "genotypes", "1", "vanraden",
               "giv 0.00", "G ASReml", "print_giv=asc", "print_geno=no genotypes.dat",
               "1", ""]
        write("tmp/calc_grm.inp", join(inp, '\n'))
        gt = open("tmp/genotypes.txt", "w")
        for i in 1:nid
            print(gt, i, ' ', join(M[:, i]), '\n')
        end
        close(gt)
    end
    G = Z'Z.*r2pq
    return p, G
end
