using Dates, HTTP
"""
    latest_beagle()
---
# Description 
- Return the latest beagle URL.
################################################################################
Note Beagle is not version controlled by git.  It seems that the names for new
versions are following the pattern like `beagle.ddMmmyy.xxx.jar`.  So I just
grab all the file names in https://faculty.washington.edu/browning/beagle/.
There must also be better ways to get the URL.
"""
function latest_beagle()
    later(a, b) = begin
        mon = Dict("Jan" => 1, "Feb" => 2, "Mar" => 3, "Apr" => 4,
                   "May" => 5, "Jun" => 6, "Jul" => 7, "Aug" => 8,
                   "Sep" => 9, "Oct" => 10, "Nov" => 11, "Dec" => 12, "mmm" => 0)
        da, db, ya, yb = parse.(Int, [a[8:9], b[8:9], a[13:14], b[13:14]])
        ma, mb = mon[a[10:12]], mon[b[10:12]]
        x = ya * 10000 + ma * 100 + da
        y = yb * 10000 + mb * 100 + db
        return x > y ? a : b
    end

    BeagleURL = "https://faculty.washington.edu/browning/beagle/"
    beagle = "beagle.00mmm00.xxx.jar"
    web = HTTP.request("GET", BeagleURL)
    txt = split(String(web.body), "\n")
    for line in txt
        r = match(r"beagle\.\d{2}\w{3}\d{2}.{5}jar", line)
        r == nothing || (beagle = later(beagle, r.match))
    end
    message("The current Beagle is $beagle.")
    return joinpath(BeagleURL, beagle)
end

"""
    latest_plink()
---
Return URL of the latest plink version 1 for Linux x86_64.
"""
function update_plink()
    message("Updating plink to the latest")
    plinkURL = "http://s3.amazonaws.com/plink1-assets"
    latest = "plink_linux_x86_64_latest.zip"
    odir = pwd()
    cd(abg_bin)
    download(joinpath(plinkURL, latest), latest)
    run(`unzip -u $latest`)
    cd(odir)
end

"""
    make_bins()
---
This is to make binaries from my C++ files.
It is supposed to run at the startup.
But this can also be run after new C++ codes were written.
"""
function make_bins()
    subtitle("C++ binaries")
    isdir(abg_bin) || mkdir(abg_bin)

    message("Compiling ABG C++ codes in `cpp` into `bin`")
    # free-space is removed, as it uses some new c++ features.
    width = 70
    for cpp in readdir(abg_cpp)
        bin = cpp[1:end-4]
        if (!isfile("$abg_bin/$bin")) || (stat("$abg_bin/$bin").mtime < stat("$abg_cpp/$bin.cpp").mtime)
            print(lpad("g++ -O3 -Wall -std=c++17 cpp/$cpp -o bin/$bin", width))
            run(`g++ -O3 -Wall -std=c++17 -o $abg_bin/$bin $abg_cpp/$cpp`)
            done()
        else
            print(lpad("bin/$bin", width))
            done("OK")
        end
    end
end

"""
    Update()
---
Every time one starts this package, the package will automatically check if `plink` and `beagle.jar` exist.
If not, the package will download the latest version of these two files.
The package will update, or download them on 17th every month anyway.
"""
function Update(force = false)
    isdir(abg_bin) || mkdir(abg_bin)
    if Dates.day(now()) == 17 || force
        rm(beagle, force=true)
        rm(plink,  force=true)
    end
    
    subtitle("Updating the binaries")

    item("Beagle related")
    beagle2vcfURL = "https://faculty.washington.edu/browning/beagle_utilities/"
    beagle2vcf = "beagle2vcf.jar"
    
    isfile(beagle) || download(latest_beagle(), beagle)
    isfile(joinpath(abg_bin, beagle2vcf)) || download(joinpath(beagle2vcfURL, beagle2vcf), joinpath(abg_bin, beagle2vcf))
    done()
    
    item("plink")
    isfile(plink) || update_plink()
    done()

    item("C++ binaries")
    make_bins()
end
