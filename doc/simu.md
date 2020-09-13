# QTL and phenotype simulation

I extracted 1000 SNP for 100 ID from a dataset I analysed.
The ID and SNP names were changed.
I also changed the chromosome and bp location of each SNP, as if they are evenly spaced on a same chromosome.

In Julia REPL:

```julia
using Manuals
qp_simulation(path=pwd(), nq=100, h2=0.8)
```

The simulation results will be saved in current path by default.
One can specify a path the results.
The number of QTL `nq` and $h^2$ can also be specified else.


