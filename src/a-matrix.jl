using DataFrames
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
    title("Calculate inversed A matrix with pedigree $ped")
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
    title("Calculate inbreeding coefficients with pedigreee $ped")
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
    title("Calculate A matrix with pedigree $ped")
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

"""
    sort_pedigree(pedigree::AbstractString, miss::AbstractString="0")
---
Given a pedigree file, this function will recode the pedigree to consective integers start from `1`.
Input file should have at least 3 columns.  The first 3 columns are `ID`, `sire` and `dam`.
ID in sire and dam columns are coded first.
Then the ID in ID columns are coded if their parents were already coded.
This procedure stops when all ID are coded, or a `loop`, i.e., an offspring is also its ancestor.
In the latter case, an error is reported.

By default, string `0` is for missing.  One can specify other string according to your data.

A 3-column DataFrame is returnd.  Row numbers are the recoded integers.  The 3rd column are their original
names.  The 1st and 2nd columns are their sires and dams.
"""
function sort_pedigree(pedigree::AbstractString, miss::AbstractString="0")
    title("Sort and recode pedigre $pedigree")
    oped = DataFrame(ID = String[], Sire = String[], Dam = String[]) # The pedigree to be sorted
    for line in eachline(pedigree)
        id, pa, ma = split(line)[[1,2,3]] # manual input meaning pedigree can have more columns
        push!(oped, (id, pa, ma))
    end

    coded = Dict()
    push!(coded, miss => 0)     # Name::String and its coded integer

    iid = 1                                                      # The code starts from 1
    sped = DataFrame(Sire = Int[], Dam = Int[], Name = String[]) # The sorted pedigree

    # ID only in the Sire and Dam column
    Parents = Set{String}()
    for id in [oped.Sire; oped.Dam]
        id == miss || push!(Parents, id)
    end
    for id in setdiff(Parents, oped.ID)
        push!(sped, (0, 0, id))
        push!(coded, id => iid) # ID in parent columns only are coded first
        iid += 1
    end

    left = prev = size(oped, 1)
    while left > 1              # got through the pedigree
        for (id, pa, ma) in eachrow(oped)
            if haskey(coded, pa) && haskey(coded, ma)
                push!(sped, (coded[pa], coded[ma], id))
                push!(coded, id => iid)
                iid += 1
                left -= 1
            end
        end
        if left == prev         # left cant't be reduced but >0
            warning("!!! Loop in the pedigree")
            break
        end
        prev = left
    end
    return sped
end
