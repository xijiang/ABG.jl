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
    # free-space is removed, as it uses some new c++ features.
    bins = ["amat", "raw2gt", "freq", "read-block"]
    width = 70
    for bin in bins
        if (!isfile("$abgBin/$bin")) || (stat("$abgBin/$bin").mtime < stat("$abgCpp/$bin.cpp").mtime)
            print(lpad("g++ -O2 -Wall -std=c++17 cpp/$bin.cpp -o bin/$bin", width))
            run(`g++ -O2 -Wall -std=c++17 -o $abgBin/$bin $abgCpp/$bin.cpp`)
            done()
        else
            print(lpad("bin/$bin", width))
            done("OK")
        end
    end
    
    # blas = ["bigm"]
    # for bin in blas
    #     if (!isfile("$abgBin/$bin")) || (stat("$abgBin/$bin").mtime < stat("$abgCpp/$bin.cpp").mtime)
    #         print(lpad("g++ -O2 -Wall -lopenblas -lpthread -std=c++17 cpp/$bin.cpp -o bin/$bin", width))
    #         run(`g++ -O2 -Wall -std=c++17 -lopenblas -lpthread -o $abgBin/$bin $abgCpp/$bin.cpp`)
    #         done()
    #     else
    #         print(lpad("bin/$bin", width))
    #         done("OK")
    #     end
    # end
end
