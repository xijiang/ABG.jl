"""
    pedsort(src, unknown, out)
---

# Description
Given a input stream `src` contains columns:
1. id
2. sire
3. dam
4. other information

And a `Set` of `unknown` strings, this program will sort the pedigree, such that `Sire` and `Dam` appear in
the `ID` column first.  In other words, an offspring appears after its parents in the sorted pedigree.

Results will be written to `out`, which has columns:
1. recoded `sire` name
2. recoded `dam` name
3. original name of current `ID`, which is the current row number (starts from 1)
4. `other information` in the orginal stream.

All former `unknown`s are coded as `0`.
"""
function pedsort(src, unknown::Set{String}, out)
    title("Sort pedigree in $src")
    #for line
end
