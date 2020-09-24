"""
    DMU()
---
# Description
Prepare data and driver files for DMU(https://dmu.ghpc.au.dk/).
One can prepare a sample dataset with `ABG.simu_qp(h2=0.8)`.
Since there is no pedigree info in toy dataset, that part is ignored here.

# Parameters

Below are for DMU1, a must-run procedure of DMU.

## Analyse
- `task` = 1, variance components estimation using REML
- `method` = 31, AI (average information), updates ⊄ par space && AI + EM
- `scaling` = 0, no data scaling before computation
- `test_prt` = 0, standard

## Data
- `FMT` = ASCII
- (#`int`, #`real`, `miss`) = (1, 1, -99)
- `fn` = project.dat, file name

## VARIABLE
- ID phntp

## MODEL

## VAR_STR
- `r_factor`: 1, structure number
- `type`: GREL, 
- `FMT`: ASCII
- `fn`: toy.giv

$VAR_STR 1 GREL ASCII /home/theom/data/ulf/4k_qtl/10000genots/chr1/mcmc2/melb/dmu_/g.giv  
"""
function DMU(;
             dir=abg_rst,
             project="toy",
             desc=["Analysis on toy data"])
    subtitle("Prepare $project data and driver files for DMU")
    item("The driver file")
    driver = joinpath(dir, "toy.drv")
    
    drvctt = ["\$COMMENT",
              join(desc, '\n'),
              "",
              "\$ANALYSE 1 31 0 0", # task=1, var
              "",
              "\$DATA ASCII (2, 2, -99) $project.dat",
              "",
              "\$VARIABLE",
              "ID phntp",
              "",
              "\$VAR_STR",
              "",
              ]
    write(driver, join(drvctt, '\n'))
    open(joinpath(dir, "toy.dat"), "w") do io
        id = 1
        for line in eachline(joinpath(dir, "toy.ph"))
            println(io, id, ' ', line)
            id += 1
        end
    end
    done()
end

#= sample driver for DMU1
$COMMENT
Dynamisk hy-definisjon
Analyse av kg melk, 1-3 laktasjon, gjentak-dyremodell random HY gen grup

$ANALYSE 1 31 0 0

$DATA ASCII (2,2,-99) fat.dat

$VARIABLE
# 1 2
id  line
#1 2   
yld wg 


$MODEL
1
0
1 2 2 2 1
1 1
0
0

$VAR_STR 1 GREL ASCII /home/theom/data/ulf/4k_qtl/10000genots/chr1/mcmc2/melb/dmu_/g.giv  

$PRIOR 
1 1 1 139666
2 1 1 224499


$DMUAI
10
1.0d-7
1.0d-6
1
0
=#
