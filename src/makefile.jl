"""
    make()
---
This is to make binaries from my C++ files.
It is supposed to run at the startup.
But thic can also be run after new C++ codes were written.
"""
function make()
    title("C++ binaries")
    isdir(abgBin) || mkdir(abgBin)

    message("Compiling ABG C++ codes in `cpp` into `bin`")
    bins = ["dnt",
            "inbreeding-coefficient"]
    for bin in bins
        if (!isfile("$abgBin/$bin")) || (stat("$abgBin/$bin").mtime < stat("$abgCpp/$bin.cpp").mtime)
            print(lpad("g++ -O2 -Wall -std=c++17 -o bin/$bin cpp/$bin.cpp", 50))
            run(`g++ -O2 -Wall -std=c++17 -o $abgBin/$bin $abgCpp/$bin.cpp`)
            done()
        else
            print(lpad("bin/$bin", 50))
            done("OK")
        end
    end
end
