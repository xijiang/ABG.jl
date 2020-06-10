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

"""
    big_G(file, target, memory_size, nthread)
---
# Introduciton
This procedure is to calculate G-matrix with HD genotypes, where N_id << N_locus.

# Data preparation
The genotypes can be prepared from plink.bed.  A small c++ program `raw2gt` is also
available with this package.

```bash
plink --cow --bfile data.bed --recode A --out target # produce target.raw
cat target.raw | raw2gt > raw.gt
```

# Example
`ABG.big_G("raw.gt", target, 10, 12)`

- **`raw.gt`** is the file prepared above
- **`target`** e.g., `target.G` is the calculation results
- **`10`**, is the amount of memory to be used in `gigabytes`. This amount is to hold two genotype blocks.
  The program won't use more than 4/5 of available machine memory.
- **`12`**: Number of threads to be used. 

This program uses vanRaden method I.
"""
function big_M(file::AbstractString, target::AbstractString, msize::Float64, nthread::Int64)
    title("Calculate G matrix with file $file, and vanRaden method I")
    tmem = round(Sys.total_memory()/2^30; digits = 2)
    tthread = Sys.CPU_THREADS
    fsize = stat(file).size
    nlc = length(split(readline(file)))
    nid = Int(fsize/nlc/2)
    nline = Int(floor(msize / 2 * 1024^3 / 8 / nlc))
    nblk = Int(ceil(nid/nline))
    free = split(read(pipeline(`$abgBin/free-space .`), String))[2]
    free = round(parse(Float64, free)/1024^3; digits=2)
    msg = lpad("System total memory: ", 36) * "$tmem GiB\n" *
        lpad("Memory to be used: ", 36) * "$msize GiB\n" *
        lpad("Total system threads: ", 36) * "$tthread\n" *
        lpad("Number of threads to be used: ", 36) * "$nthread\n" *
        lpad("File size: ", 36) * "$fsize bytes\n" *
        lpad("Number of loci: ", 36) * lpad("$nlc\n", 8) *
        lpad("Number of ID: ", 36) * lpad("$nid\n", 8) *
        lpad("Number of ID to deal a time: ", 36) * lpad("$nline\n", 8) *
        lpad("Number of blocks: ", 36) * lpad("$nblk\n", 8) *
        lpad("Free space in current path: ", 36) * "$free GiB\n" *
        lpad("Tmp file usage: ", 36)
    message(msg)
    i = 0
    for line in eachline(file)
        t = split(line)
        i += 1
        if i == 1000
            break
        end
    end
end
