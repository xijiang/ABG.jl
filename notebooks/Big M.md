# G matrix calculation with HD genotypes
## Problem description
Calculation of genomic relationship matrices (GRM) is a commonplace in breeding practices.
When the population is genotyped with high density (HD) measures, the calculation of GRM became difficult as it is even difficult to store the entire genotype results into memory.
For example, the memory needed for a genotyped population of 60000 ID, and the chip density 600000 SNP is $\frac {60000\times 600000\times 8}{1024^3}=268.22$ GiB.
Few computers have such amount of memory.
It is then convenient to divide the genotype matrix into several blocks of ID groups, such that a computer can hold two of such blocks a time for the block-wised calculation of GRM.

Another problem is for the general matrix multiplication (GEMM).
Although the operation is rudimentary in linear algebra, GEMM has been gone though intensive study.
Doing GEMM simply by definition is unacceptable for such situation.
For example, if we calculate GRM as below:

```julia
for i in 1:nid
    for j in 1:i
    	for k in 1:nlc
	    	GRM[i,j] += M[i,k] * M[j,k]
		end
    end
end
```

It may take months for the data of above dimension on most servers.
Using mature libraries, e.g., OpenBLAS, may shorten the calculation thousands of times.
It is however difficult to tune the parameters the library required.
For example in C OpenBLAS, it is written blow:
```C
cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, m, nlc, alpha, &mi[0], nlc, &mi[0], nlc, beta, &gblk[0], m);
```
for a diagonal block.

In the language of Julia, this is much simpler, e.g,
```julia
using LinearAlgebra 
BLAS.set_num_threads(12)
m'm  # GEMM as simple as this
```
It is also much easier to specify how many CPU threads for the calculation.

So the object of this function is to use Julia for a block-wised GRM calculation.

## Usage
If you are already using Julia, then type ']' in Julia REPL to enter the package management environment:
```julia
add https://github.com/xijiang/ABG.jl
```
Go back to the REPL and type:
```julia
using ABG # for Algorithms in Breeding and Genetics
AGB.make() # to compile some c++ codes
?ABG.G_with_big_M  # for the program documentation
ABG.G_with_big_M("gt.txt", "freq.txt", 20., 12, "out.G")
```
Above will calculate a GRM, 'out.G', with genotypes in 'gt.txt', and allele frequencies specified in 'freq.txt', using vanRaden method I.
Other methods can be integerated upon request.
The '20.' is the memory (in GiB) required for the calculation.
Note, this must be a float.
If you just type '20', there will be an error prompted.
The parameter '12' is the number of threads used for the calculation.

## Data preparation
It is convenient to prepare the genotypes and frequencies from plink files.
Suppose the genotypes are stored in plink.{bed,bim,fam}:
```bash
plink --species --bfile plink --recode A --out gt
cat gt.raw | pathto/raw2gt freq.txt > gt.txt
rm gt.raw # to save disk space
```
Program `raw2gt` came with this package.
It removes the first 6 columns from gt.raw, and all spaces except the newline character.
It also calculate the allele frequencies to freq.txt as specified above.
In the end, the 012 genotypes for each ID are in there own line, with no spaces.

### Some tweaks
#### Other frequencies
You can also specify other frequencies, e.g.,:
```bash
yes 0.5 | head -nid >freq.txt
```
to use frequency 0.5 for all loci.

Note also that this funciton does not modify the diagonals.
In the case to make the GRM positive definite, you can
#### Diagonals
```julia
using LinearAlgebra
G += Diagonal(ones(nid).*1e-6) # for example
```
#### Output formats
Currently the results are still in binary blocks. Comments or output format suggestions are welcome.

### Problem
- Validation
Bug report are heartfully appreciated.

- Frequency
Package plink can also calculate frequcies taking the advantage of binary storage. 
It however only calculate the frequencies in the founders.

## Performance
With 12 threads of AMD 3900x, ~20GiB memory and a 500GB NVMe SSD, this funciton finish the GRM of 60745 ID and 62219 loci n 2 hours.