"""
    G_with_big_M(file, target, memory_size, nthread)
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
function G_with_big_M(gt::AbstractString, frqs::AbstractString, msize::Float64, nthread::Int64)
    title("Calculate G matrix with file $gt, and vanRaden method I")
    item("Prepare the system")
    tmem = round(Sys.total_memory()/2^30; digits = 2)
    tthread = Sys.CPU_THREADS
    asked = nthread
    if nthread > tthread
        warning("No enough physical threads for the asked\nNumber of thread was set to 1")
        nthread = 1
    end
    BLAS.set_num_threads(nthread)
    fsize = stat(gt).size
    nlc = length(readline(gt))
    lsz = nlc + 1
    nid = Int(fsize/(nlc+1))
    if nid * (nlc+1) != fsize
        warning("Not a square file")
        return
    end
    frq = Float64[]
    for f in eachline(frqs)
        push!(frq, parse(Float64, f))
    end
    if length(frq) ≠ nlc
        warning("Number of loci in $frqs doesn't match")
        return
    end
    twop = 2 .* frq

    nbid = Int(floor(msize / 2 * 1024^3 / 8 / nlc * .7)) # .7 is arbitrary to avoid GC problem
    nblk = Int(ceil(nid/nbid))
    free = split(read(pipeline(`$abgBin/free-space .`), String))[2]
    free = round(parse(Float64, free)/1024^3; digits=2)
    disk = round(nid * nid * 8. /1024^3 * 2; digits=2)
    
    if disk >= free
        message("Not enough disk space")
        return
    end
    done()

    blks = [range(   1, step=nbid, length=nblk) range(nbid, step=nbid, length=nblk)]
    blks[end] = nid             # define blocks
    
    item("Summary of your system and parameters")
    msg = lpad("System total memory: ", 36) * "$tmem GiB\n" *
        lpad("Memory to be used: ", 36) * "$msize GiB\n" *
        lpad("Total system threads: ", 36) * "$tthread\n" *
        lpad("Number of threads asked: ", 36) * "$asked\n" *
        lpad("Number of threads to be used: ", 36) * "$nthread\n" *
        lpad("File size: ", 36) * "$fsize bytes\n" *
        lpad("Number of loci: ", 36) * lpad("$nlc\n", 8) *
        lpad("Number of ID: ", 36) * lpad("$nid\n", 8) *
        lpad("Number of ID to deal a time: ", 36) * lpad("$nbid\n", 8) *
        lpad("Number of blocks: ", 36) * lpad("$nblk\n", 8) *
        lpad("Free space in current path: ", 36) * lpad("$free GiB\n", 12) *
        lpad("Result files usage: ", 36) * lpad("$disk GiB\n", 12)
    message(msg)
    write("tmp/summary.txt", msg)

    item("Calculate the blocks")
    fgt = open(gt, "r")
    function readGt(x, y)
        t = Float64[]
        seek(fgt, (x-1)*lsz)
        for k in x:y
            line = readline(fgt)
            for c in line
                if c == '0'
                    push!(t, 0.)
                elseif c == '1'
                    push!(t, 1.)
                else
                    push!(t, 2.)
                end
            end
        end
        reshape(t, :, y-x+1)
    end

    for i in 1:nblk
        x, y = blks[i, :]
        print('\r', lpad("Block $i x $i of $nblk", 50))
        mi = readGt(x, y) .- twop
        write("tmp/blk-$i-$i.bin", mi'mi)
        for j in i+1:nblk
            a, b = blks[j, :]
            print('\r', lpad("Block $i x $j of $nblk", 50))
            mj = readGt(a, b) .- twop
            write("tmp/blk-$i-$j.bin", mi'mj)
        end
    end
end
