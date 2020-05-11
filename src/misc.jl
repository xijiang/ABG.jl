
"""
    empty_dir(dir::AbstractString = "tmp")
---
If exists `dir`, empty everything in it.  Or just return.
"""
function empty_dir(dir::AbstractString = "tmp")
    isdir(dir) || return
    for f in readdir(dir)
        rm(joinpath(dir, f), force=true, recursive=true)
    end
end
