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

Each execution mode computes and generates a combination of files.

#### Depth of coverage characterization of all genomic regions
For each region provides the mean depth of coverage, the GC content and the mean depth of coverage of the subregion (user specified, default 0.5 fraction) that maximizes the coverage peak signal (`rcS` and correspodning genomic coordinates `fromS` and `toS`), to account for the reduced coverage depth due to incomplete match of reads to the captured regions.

```
chr	from	to	fromS	toS	rc	rcS	gc
1	26823507	26823727	26823564	26823673	412.30	497.55	0.60
1	26825032	26825177	26825062	26825133	711.43	754.71	0.42
1	26826977	26827132	26827000	26827076	180.24	198.53	0.72
1	26828102	26828248	26828136	26828208	492.88	529.68	0.58
1	26828987	26829053	26829011	26829043	595.55	610.88	0.42
1	26836987	26837058	26836994	26837028	161.72	164.60	0.59
...
```

## Visual reports

## Licence

## Credit
