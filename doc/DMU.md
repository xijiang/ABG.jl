# Notes for the DMU software

## Download and installation
### the website

The software can be downloaded from https://dmu.ghpc.au.dk/DMU/.  The following notes
are just for Linux version.

### Install the Linux version

Download the tar ball from https://dmu.ghpc.au.dk/DMU/Linux/Current/.  E.g.,

```bash
file="dmuv6-R5.3-EM64T-build-2019-02-08.tar.gz"
curl -O https://dmu.ghpc.au.dk/DMU/Linux/Current/$file
tar xvf $file
mkdir -p ~/.local/bin
cp dmuv6/R5.3-EM64T/bin/* ~/.local/bin
```

## Example for often used functions
### 