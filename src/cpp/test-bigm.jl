"""
This function
- generate a random 012 genotype set
- write them to a file
- calculate t't
- call bigm.cpp to do the dame calculation
- compare results to make sure the results are right
"""
function test_bigm()
    #= create a test set
    nid, nlc = 1000, 20000
    ndw = 2 * nid
    t = rand(0:2, nlc, nid)
    open("tmp/genotypes", "w") do io
        for i in 1:nid
            println(io, join(t[:, i]))
        end
    end
    frq = sum(t, dims=2) ./ ndw
    open("tmp/freq.txt", "w") do io
        for x in frq
            println(io, x)
        end
    end
    =#
end
