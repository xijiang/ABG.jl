"""
    A_inverse(ped::AbstractString)
---
################################################################################
Given a 2-column pedigree file `ped`, this function returns a SparseArray of `A`
matrix inversed.

# Pedigree file format
- Each line has two `Int` numbers for `sire` and `dam`, respectively and in
  that order.
- Number `0` is for missing parents.
- Line numbers start from `1` are deemed as ID number.

See `dat/henderson1976.ped` for example.
"""
function A_inverse(ped::AbstractString)
    Ai = Any
    open(pipeline(`cat $ped`, `$abgBin/amat D`), "r") do io
        nid = parse(Int, readline(io))
        D = Float64[]
        for _ in 1:nid
            d = parse(Float64, readline(io))
            push!(D, d)
        end
        Ti = Int[]
        Tj = Int[]
        Td = Float64[]
        for line in eachline(io)
            i, j, d = split(line)
            push!(Ti, parse(Int, i))
            push!(Tj, parse(Int, j))
            push!(Td, parse(Float64, d))
        end
        T = sparse(Ti, Tj, Td)
        Ai = T'Diagonal(D)T
    end
    return Ai
end

"""
    inbreeding_coefficient(ped::AbstractString)
---
Given a 2-column pedigree file `ped`, this program calculate inbreeding
coefficient of everybody in the pedigree.
"""
function inbreeding_coefficient(ped::AbstractString)
    ic = Float64[]
    open(pipeline(`cat $ped`, `$abgBin/amat F`), "r") do io
        for line in eachline(io)
            t = parse(Float64, line)
            push!(ic, t)
        end
    end
    return ic
end

"""
    A_matrix(ped::AbstractString, list::AbstractString="")
---
Given a 2-column pedigree file `ped`, this function calculate an **`A`** matrix.
IF a subset of ID `list` of the predigree is given, this gives a partial
**`A`** matrix. Inverse of this partial **`A`** can be done within `Julia`.
"""
function A_matrix(ped::AbstractString, list::AbstractString="")
    ID = Int[]
    A = Float64[]
    println(ped)
    open(pipeline(`cat $ped`, `$abgBin/amat A $list`), "r") do io
        nid = parse(Int, readline(io))
        for _ in 1:nid
            id = parse(Int, readline(io))
            push!(ID, id)
        end
        for line in eachline(io)
            for elm in split(line)
                push!(A, parse(Float64, elm))
            end
        end
        reshape(A, (nid, nid))
    end
    return ID, A
end
