## The `plink` options I often used

- Don't forget `--species=cow`, for example, for animal breeders.

### Convert `ped + map` $\to$ `bed`

```bash
--recode --make-bed --ped $ped --map $map --out $out
```

### The I/O `bed` dataset

```bash
--bfile $src			# In dataset
--make-bed --out $out		# Out dataset
```

- I will ignore these below.

### Extract a subset of SNP
```bash
--extract $snp			# snp is the file contains SNP names.
```

### Missing allele statistics

```bash
--chr $chr --missing		# chr can be, e.g., `1-29`.
```

### MAF statistics

```bash
--nonfounders --freq		# --nonfounders is to include everybody.
```

### Hardy-Weinberg test

```bash
--chr $chr --hardy
```

### Merge several `plink` datasets

```bash
--merge-list $list		# list is the file name contains plink dataset name stems
```

### Convert to VCF

```bash
--recode vcf-iid bgz		# vcf-iid is to ignore family info;  bgz is to save to `.gz`
```

### VCF to plink

```bash
--vcf $src
```

### Convert to `map+ped`

```bash
--recode
```

### QC

```bash
--geno $geno			# 0.1
--maf $maf			# 0.01
--hwe $hwe			# 1e-4
--mind $mind			# 0.1
```

### ID subset

```bash
--keep $list			# file list has two columns: family ID
```

### Convert to 012 genotypes

```bash
--recode A
```
