"""
    readIDnGt(file::AbstractString)
---
Read ID and genotypes from `file`.  The first column in `file` are ID names.  The rest 
columns are genotypes for each ID on every locus.

No error check.
"""
function readIDnGt(file::AbstractString)
    title("ID names and genotypes")
    ID = String[]
    Z  = Float64[]
    for line in eachline(file)
        v = split(line)
        push!(ID, v[1])
        append!(Z, parse.(Float64, v[2:end]))
    end
    return ID, reshape(Z, :, length(ID))
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
function vanRaden(Z::Array{Float64, 2})
    title("Calculate G with vanRaden method I")
    item("Allele frequencies")
    twop = mean(Z, dims = 2)
    sum2pq = twop' * (1 .- .5 .* twop)
    Z .-= twop
    G = Z'Z ./ sum2pq
end

"""
    vanRaden2(W::Array{Float64, 2})
---
Calculate **`G`** with van Raden method II.
"""
function vanRaden2(W::Array{Float64, 2})
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

    
