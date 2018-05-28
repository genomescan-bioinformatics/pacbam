# PaCBAM

PaCBAM is a C command line tool for the complete characterization of genomic regions and single nucleotide positions from next-generation sequencing data.  
PaCBAM implements a fast and scalable multi-core computational engine, generates exhaustive output files for downstream analysis, introduces an innovative on-the-fly read duplicates filtering strategy and provides comprehensive visual reports. 

## Installation

To install PaCBAM clonet the repository and compile the C source code.

```bash
git clone https://CibioBCG@bitbucket.org/CibioBCG/pacbam.git  
cd pacbam
make -f Makefile.linux
```

Use instead Makefile.osx and Makefile.mingw to compile PaCBAM for, respectively, OSX and Windows systems.  
Samtools library `libbam.a` has been generated for GNU/Linux, Windows and MacOSX systems.  
For compilation on Windows we have added also `libz.a` library, while compilation on Linux/MacOSX requires the installation of the development `zlib` package.  
Libraries can be found in `./lib` directory.  
Windows libraries have been generated using MinGW.  
If libraries are not working we suggest to download them again.  

## Usage

## Output files

## Visual reports

## Licence

## Credit
