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

Each execution mode computes and generates a combination of the followinf files.

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
1	26843222	26843294	26843244	26843279	423.22	436.81	0.40
1	26844137	26844233	26844163	26844210	183.45	191.83	0.47
1	26845287	26845429	26845323	26845393	952.87	1026.44	0.50
...
```

#### Single-base resolution pileup

For each genomic position in the target provides the read depth of the 4 possible bases A, C, G and T, the total depth of coverage, the variants allelic fraction (VAF), the strand bias information for each base, the unique identifier (e.g. dbsnp id) if available.

```
chr	pos	ref	A	C	G	T	af	cov
1	26823507	C	0	186	0	0	0.000000	186
1	26823508	A	167	0	0	0	0.000000	167
1	26823509	G	0	0	166	0	0.000000	166
1	26823510	C	0	180	0	0	0.000000	180
1	26823511	T	0	4	0	177	0.022099	181
1	26823512	G	0	0	187	0	0.000000	187
1	26823513	C	0	195	0	0	0.000000	195
1	26823514	A	196	0	0	0	0.000000	196
1	26823515	G	0	0	202	0	0.000000	202
```

## Visual reports

## Licence

## Credit
