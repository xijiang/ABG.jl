# ABG.jl
Algorithms for Breeding and Genetics

I don't export any function in this package. But I have defined some shorthand
examples that can be included in your own project, e.g., `v1.1.jl`.

This package will create some directory and files relative to the path you
imported this package.

# Version history
## v0.2
- $\mathbf{A}$ matrix related.

## v0.1
- Styled printing
- Frequently used `plink` commands
  - You **must** specify where to find `plink` when you use these commands
  - For example: `plink = "bin/plink"`