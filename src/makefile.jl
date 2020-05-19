"""
    make()
---
This is to make binaries from my C++ files.
It is supposed to run at the startup.
But thic can also be run after new C++ codes were written.
"""
function make()
    title("C++ binaries")
    pkg = begin
        t = splitdir(pathof(ABG))[1]
        l = findlast('/', t)
        SubString(t, 1:l)
    end
    src = joinpath(pkg, "cpp")
    bin = joinpath(pkg, "bin")
    isdir("bin") || mkdir("bin")

    message("Compiling ABG C++ codes in $src into $bin")
    #bins = []
    # for bin in bins
    #     if (!isfile("bin/$bin")) || (stat("bin/$bin").mtime < stat("c/$bin.cpp").mtime)
    #         print("g++ -O2 -Wall -o bin/$bin src/$bin.cpp")
    #         run(`g++ -O2 -Wall -o bin/$bin src/$bin.cpp`)
    #         done()
    #     else
    #         print("bin/$bin")
    #         done("OK")
    #     end
    # end
end
