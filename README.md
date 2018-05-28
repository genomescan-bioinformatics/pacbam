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

Use instead Makefile.macosx and Makefile.mingw to compile PaCBAM on, respectively, OSX and Windows systems.  
Samtools library `libbam.a` has been generated for GNU/Linux, Windows and MacOSX systems.  
For compilation on Windows we have added also `libz.a` library, while compilation on Linux/MacOSX requires the installation of the development `zlib` package.  
Libraries can be found in `./lib` directory.  
Windows libraries have been generated using MinGW.  
If libraries are not working we suggest to download/recompile them again.

## Usage

Running PaCBAM executable will list usage options. 

```
Usage: 
 ./pacbam bam=string bed=string vcf=string fasta=string [mode=int] [threads=int] [mbq=int] [mrq=int] [mdc=int] [out=string] [duptab=string] [regionperc=float] [strandbias]

bam=string 
 NGS data file in BAM format 
bed=string 
 List of target captured regions in BED format 
vcf=string 
 List of SNP positions in VCF format 
fasta=string 
 Reference genome FASTA format file 
mode=string 
 Execution mode [0=RC+SNPs+SNVs|1=RC+SNPs+SNVs+PILEUP(not including SNPs)|2=SNPs|3=RC|4=PILEUP]
 (default 0)
duptab=string 
 On-the-fly duplicates filtering lookup table
threads=int 
 Number of threads used (if available) for the pileup computation
 (default 1)
regionperc=float 
 Fraction of the captured region to consider for maximum peak signal characterization
 (default 0.5)
mbq=int 
 Min base quality
 (default 20)
mrq=int 
 Min read quality
 (default 1)
mdc=int 
 Min depth of coverage that a position should have to be considered in the output
 (default 0)
strandbias 
 Print strand bias count information
 (default 1)
out=string 
 Path of output directory (default is the current directory)
```

## Output files

## Visual reports

## Licence

## Credit
