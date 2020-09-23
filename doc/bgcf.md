# User manual for bgcf
This program was written by Theo.  I obtained from him directly.

## Driver file
This file should be called `bgcf.inp`, and have the following contents:

| Item | Comment |
| -- | -- |
| &inputvalues | |
| Nanim=35688 | !nanim=number of animals in data |
| Nchr=1 | !nchr=number of chroms (limited to 1 for now) |
| Nsnptotal=120160 | !nsnptotal |
| bimfil=chr29_noxbred.bim | !PLINK bim file (if nchr>1 => chromosome nr will be appended to filename) |
| bedfil=chr29_noxbred.bed | !PLINK bed file (if nchr>1 => chromosome nr will be appended to filename) |
| phenfil=hbc.fatpr | !phenfilename of phenotype file containing nanim rows and 3 columns: Y1  weight fixed_effect_class (weight=1 => equal weights) |
| PI_value=.002 | !PI-value for BayesC analysis (between 0 and 1) (negative PI =>  starting value for estimation of PI) |
| Ve=1.31e-5 | !R=error covariance matrix across traits (full-matrix is needed; only this last line may contain comment; same other matrices)  |
| Vsnp=2.26e-7 | !Vsnp matrix; SNP effect ~ N(0,Vsnp) 1% of Vg |
| Vu=2.26e-5 | !G = polygenic covariance matrix (all 0's indicates no polygenic variance) 70% of Vg |
| Ncycle=1000 | !ncycle=number of gibbs cycles |
| Nburnin=100 | !nburnin=number of burnin cycles |
| Nwork=1 | !number of parallel chains == number of threads to be used (dont use more threads then available on computer) |
| Ithin=0 | !thinning output (print every ithin sample) (ithin=0 => no output of samples (only summary statistics)) |
| Ispeedup=0 | !Ispeedup; reduce computations by factor Ispeedup (Ispeedup=0 => no reduction) |
| / | |

## Notes:
- The file can contain only the first column, starts from `&`, and stops at `/`. 
- Parameters like $V_e$, and $V_u$ can be estimated with `DMU`.  See `DMU.md` came with this package on its commonly used functions.
- Set $\text{PI}_{\text{value}$ and $V_{\text{snp}}=0$ for a no QTL model.

## Prepare data and the driver file for bgcf

```julia
# ]
# add https://github.com/xijiang/ABG.jl
# ^h
using ABG
ABG.bgcf()
```

The result data will be in `$PWD/sandbox`.
One can change this path by calling:

```julia
ABG.set_rst("your_desired_path")
```
