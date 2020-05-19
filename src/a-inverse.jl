"""
    A_inv(ped::AbstractString)
---
Given a 2-column pedigree file `ped`, this function returns a SparseArray of `A` matrix inversed.

# Pedigree file format
- Each line has two `Int` numbers for `sire` and `dam`, respectively and in that order.
- Number `0` is for missing parents.
- Line numbers start from `1` are deemed as ID number.

See `dat/henderson1976.ped` for example.
"""
function A_inv(ped::AbstractString)
    open(ped, "r") do io
        for line in eachline(io)
            println(line)
        end
    end
end

# D   = readdlm("D.vec")[:,1]
# dat = readdlm("T.mat")
# T   = sparse(Int.(dat[:,1]), Int.(dat[:,2]), dat[:,3])
# dat = nothing
# 
# Ai  = T'Diagonal(D)T
# 
# # An efficient way for data I/O
# # write the inverse to a file, say inverse.bin
# file = open("inverse.bin", "w")
# serialize(file, Ai)
# close(file)
# 
# # read them back
# Ai = deserialize(open("inverse.bin"))
