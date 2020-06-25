"""
    G_with_big_M(gt_file, frq, memory_size, nthread, out)
---
# Introduciton
This procedure is to calculate G-matrix with HD genotypes, where N_id << N_locus.

# Data preparation
The genotypes are prepared from plink.bed.  A small c++ program `raw2gt` was
written together with this package to convert `target.raw`.
The genotype file should have (`012`) genotypes one ID a line.
The newline should of *nix convention.

```bash
plink --cow --bfile data.bed --recode A --out target # produce target.raw
cat target.raw | raw2gt raw.frq > raw.gt
```

Above can also be finished with function `plink_2_gt` of this package, e.g.,
```julia
ABG.plink_2_gt("path-to-plink-file-stem")
```
The results will be written to `plk.gt` and `plk.frq` in path `tmp` in the current directory.

# Example
`ABG.big_G("raw.gt", "raw.frq", 10., 12, "raw.G")`

- **`raw.gt`** is the file prepared above
- **`raw.G`** is the calculation results.  3-column: `raw column element-value` of the `G`-matrix.
- **`10`**, is the amount of memory to be used in `GiB`.  This amount is to hold two+ genotype blocks.
  The program won't use more than 4/5 of available machine memory.
- **`12`**: Number of threads to be used.  If you asked more than physical threads, it is set to `1`.

This program uses vanRaden method I.
"""
function G_with_big_M(gt::AbstractString,
                      frqs::AbstractString,
                      msize::Float64,
                      nthread::Int64,
                      out::AbstractString)
    title("Calculate G matrix with file $gt, and vanRaden method I")
    tmp = joinpath(workdir, "tmp")
    isdir(tmp) || mkdir(tmp)
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
    fszg = round(fsize/1024^3; digits=2) # file size in GiB
    nlc = length(readline(gt))
    lsz = nlc + 1               # assuming *nix file newlines.
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
    twop = frq.*2.
    rs2pq = begin                # !!! multiply is faster than divide
        q = 1 .- frq
        s = q'frq
        .5/s
    end
    
    nbid = Int(floor(msize / 2 * 1024^3 / 8 / nlc * .8)) # .8 is arbitrary to avoid GC problem
    nblk = Int(ceil(nid/nbid))
    #free = split(read(pipeline(`$abgBin/free-space .`), String))[2]
    free = begin
        statfs = zeros(UInt64, 11)
        ccall((:statvfs, "libm.so.6"), Cint, (Cstring, Ptr{Cvoid}), ".", statfs)
        round(convert(Int, statfs[1]*statfs[5])/1024^3; digits=2)
    end
    disk = round(nid * nid * 8. /1024^3 * 2; digits=2)

    if disk >= free
        message("Not enough disk space")
        return
    end
    done()
    blks = [range(1, step=nbid, length=nblk) range(nbid, step=nbid, length=nblk)]
    blks[end] = nid             # define blocks
    
    item("Summary of your system and parameters")
    msg = lpad("System total memory: ", 36) * "$tmem GiB\n" *
        lpad("Memory to be used: ", 36) * "$msize GiB\n" *
        lpad("Total system threads: ", 36) * "$tthread\n" *
        lpad("Number of threads asked: ", 36) * "$asked\n" *
        lpad("Number of threads to be used: ", 36) * "$nthread\n" *
        lpad("File size: ", 36) * "$fszg GiB\n" *
        lpad("Number of loci: ", 36) * lpad("$nlc\n", 8) *
        lpad("Number of ID: ", 36) * lpad("$nid\n", 8) *
        lpad("Number of ID to deal a time: ", 36) * lpad("$nbid\n", 8) *
        lpad("Number of blocks: ", 36) * lpad("$nblk\n", 8) *
        lpad("Free space in current path: ", 36) * lpad("$free GiB\n", 12) *
        lpad("Result files usage: ", 36) * lpad("$disk GiB\n", 12)
    message(msg)
    write("$tmp/summary.txt", msg)

    item("Calculate the blocks")
    fgt = open(gt, "r")
    function readGt(x, y)
        t = Array{Float64, 2}(undef, nlc, y+1-x) # column majored for easy coding later
        k = 1
        seek(fgt, (x-1)*lsz)
        for _ in x:y
            g = read(fgt, lsz)  # also read the (one) newline character
            for i in 1:nlc
                t[i, k] = Float64(g[i] - 0x30) # 0x30 == '0'
            end
            k = k+1
        end
        t .-= twop
    end

    for i in 1:nblk
        x, y = blks[i, :]
        print('\r', lpad("Block $i x $i of $nblk", 50))
        mi = readGt(x, y)
        G = mi'mi
        G .*= rs2pq
        write(joinpath(tmp, "$i-$i.blk"), G)
        for j in i+1:nblk
            a, b = blks[j, :]
            print('\r', lpad("Block $i x $j of $nblk", 50))
            mj = readGt(a, b)
            G = mj'mi
            G .*= rs2pq
            write(joinpath(tmp, "$j-$i.blk"), G)
        end
    end
    done()

    item("Merge and output to row col element-value => $out")
    # dam tedious below.
    open(out, "w") do io
        for i in 1:nblk
            x, y = blks[i, :]
            nrb = y-x+1         # number rows in the block
            for j in 1:i-1
                a, b = blks[j, :]
                ncb = b-a+1     # number columns in the block
                dblk = reinterpret(Float64, read("$tmp/$i-$j.blk", nrb*ncb*8))
                k = 1           # column majored anyway
                for col in a:b
                    for row in x:y
                        print(io, "$row $col ", dblk[k], '\n')
                        k += 1
                    end
                end
            end
            # The diagonal blocks
            dblk = reshape(reinterpret(Float64, read("$tmp/$i-$i.blk", nrb*nrb*8)), nrb, :)
            for p in x:y
                for q in x:p
                    print(io, "$p $q ", dblk[p-x+1, q-x+1], '\n')
                end
            end
        end
    end
    done()
end

"""
    plink_2_gt(plink, out="plk"; species="cow")
---
Convert plink genotypes in `plink`.bed to `012` genotypes, and then remove the first 6
columns, and spaces in the result file.

## Example
```
ABG.plink_2_gt("sample", specices="dog")
```
The results will be written to `plk.gt` and `plk.frq`,
in `tmp` of current working directory.
"""
function plink_2_gt(plink::AbstractString, out::AbstractString="plk"; species::AbstractString="cow")
    title("Convert $plink.bed to genotype file required for this package")
    message("   data: $plink.{bed,bim,fam}\nSpecies: $species")
    if Sys.which("plink") == Nothing
        warning("Command plink can't be find in default paths")
        return
    end
    isfile(joinpath(abgBin, "raw2gt")) || make()
    tmp = joinpath(workdir, "tmp")
    isdir(tmp) || mkdir(tmp)
    plink_012(plink, joinpath(tmp, out), species) # write to tmp/out.raw
    _ = run(pipeline(joinpath(tmp, "$out.raw"), `$abgBin/raw2gt $tmp/$out.frq`, "$tmp/$out.gt"))
    message("Results were written to $tmp/$out.gt")
    rm("$tmp/$out.raw", force=true)
    done()
end
