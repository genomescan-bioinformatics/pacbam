# PaCBAM

PaCBAM is a C command line tool for the complete characterization of genomic regions and single nucleotide positions from next-generation sequencing data.  
PaCBAM implements a fast and scalable multi-core computational engine, generates exhaustive output files for downstream analysis, introduces an innovative on-the-fly read duplicates filtering strategy and provides comprehensive visual reports. 

## Multi-platform binaries and Docker/Singularity containers

Binaries for Linux, macOS and Windows platforms can be found in folder `binaries`.  
Docker container can be found at ...  
Singularity container can be found at ...

## Compilation from source code

To install PaCBAM clone the repository and compile the C source code.

```bash
git clone https://CibioBCG@bitbucket.org/CibioBCG/pacbam.git  
cd pacbam
make -f Makefile.linux
```

Use instead Makefile.macos and Makefile.mingw to compile PaCBAM on, respectively, macOS and Windows systems.  
Samtools library `libbam.a` has been generated for GNU/Linux, Windows and macOS systems.  
For compilation on Windows we have added also `libz.a` library, while compilation on Linux/macOS requires the installation of the development `zlib` package.  
Libraries can be found in `./lib` directory.  
Windows libraries have been generated using MinGW.  
If libraries are not working we suggest to download/recompile them again.

## Usage
PaCBAM expects as input a sorted and indexed BAM file, a BED file with the coordinates of the genomic regions of interest (namely the target, e.g. captured regions of a WES experiment), a VCF file specifying a list of SNPs within the target and a reference genome FASTA file.  
Different running modes and filtering/computation options are available.  
Running PaCBAM executable will list all usage options. 

```
Usage: 
 ./pacbam bam=string bed=string vcf=string fasta=string [mode=int] [threads=int] [mbq=int] [mrq=int] [mdc=int] [out=string]
          [duptab=string] [regionperc=float] [strandbias]

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

## Examples

Folder `examples` contains a small example of a BAM file and correspoding target regions in BED format and a SNPs in target regions in VCF format.  
The following command executes PaCBAM with mode 1, generating 4 output files.

```bash
../pacbam bam=NGSData.bam bed=TargetRegions.bed vcf=SNPsInTargetRegions.vcf fasta=/path-to-reference-genome/human_g1k_v37.fasta mode=1 out=./
```

The reference genome to use in this example can be downloaded at

`ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz`

To activate the *on-the-fly read duplicates filtering* add to the command `duptab=FileName.txt` where file `FileName.txt` is a text tab-delimtied file with the following format:

```
 0 100 1
 100 500 2
 500 1000 3
```

This description means that the coverage of all positions with coverage in (0,100] will be normalized considering only 1 read per read duplicates group. A read duplicate group is
determined by collecting all reads at the position and grouping them considering their strand and alignment positions (both read pair position are used in paired-end sequencing experiments).  
Similar filtering strategy is used for coverage intervals (100,500] and (500,1000] where at most 2 and 3 reads per read duplicates group are used, respectively.  
Positions with coverage >1000 will use the threshold of the highest interval (i.e. 3).  
To perform complete deduplication specify a single interval with 3rd column value equal to 1. 

## Output files

Each execution mode computes and generates a combination of the following files.

#### Depth of coverage characterization of all genomic regions
For each region provides the mean depth of coverage, the GC content and the mean depth of coverage of the subregion (user specified, default 0.5 fraction) that maximizes the coverage peak signal (`rcS` and correspodning genomic coordinates `fromS` and `toS`), to account for the reduced coverage depth due to incomplete match of reads to the captured regions.

```
chr	from	to	fromS	toS	rc	rcS	gc
20	68348	68410	68348	68378	130.40	129.68	0.48
20	76643	77060	76845	77052	81.18	111.99	0.41
20	123267	123329	123293	123323	93.00	99.81	0.50
20	126053	126335	126100	126240	32.55	54.73	0.44
20	138183	138236	138210	138235	78.08	99.92	0.51
20	139412	139667	139510	139636	117.86	125.38	0.39
20	168524	168761	168524	168641	69.79	91.03	0.39
20	170213	170266	170213	170238	13.91	18.69	0.40
20	207927	207989	207958	207988	96.40	106.65	0.48
...
```

#### Single-base resolution pileup
For each genomic position in the target provides the read depth of the 4 possible bases A, C, G and T, the total depth of coverage, the variants allelic fraction (VAF), the strand bias information for each base, the unique identifier (e.g. dbsnp id) if available.

```
chr	pos	ref	A	C	G	T	af	cov
20	68348	G	0	0	129	0	0.000000	129
20	68349	C	0	130	0	0	0.000000	130
20	68350	C	0	130	0	0	0.000000	130
20	68352	T	0	0	0	130	0.000000	130
20	68353	G	0	0	130	0	0.000000	130
20	68354	A	130	0	0	0	0.000000	130
20	68355	A	130	0	0	0	0.000000	130
20	68356	T	0	0	0	130	0.000000	130
20	68357	A	130	0	0	0	0.000000	130
...
```

#### Single-nucleotide variants (SNVs) pileup
Provides pileup information only for position with positive VAF, computed using the alternative base with highest read depth (if any).

```
chr	pos	ref	alt	A	C	G	T	af	cov
20	76953	G	A	1	0	99	0	0.010000	100
20	126263	C	T	0	26	0	1	0.037037	27
20	139484	A	G	156	0	1	0	0.006369	157
20	139557	A	G	99	0	1	0	0.010000	100
20	139570	C	A	1	171	0	0	0.005814	172
20	139622	C	A	1	135	0	0	0.007353	136
20	168728	T	A	56	0	0	0	1.000000	56
20	209986	A	T	227	0	0	2	0.008734	229
20	210097	C	T	0	82	0	1	0.012048	83
...
```

when `strandbias` option is used, the output format is the following:

```
chr	pos	ref	alt	A	C	G	T	af	cov	Ars	Crs	Grs	Trs
20	76953	G	A	1	0	99	0	0.010000	100	1	0	80	0
20	126263	C	T	0	26	0	1	0.037037	27	0	0	0	0
20	139484	A	G	156	0	1	0	0.006369	157	111	0	1	0
20	139557	A	G	99	0	1	0	0.010000	100	39	0	0	0
20	139570	C	A	1	171	0	0	0.005814	172	0	91	0	0
20	139622	C	A	1	135	0	0	0.007353	136	0	67	0	0
20	168728	T	A	56	0	0	0	1.000000	56	19	0	0	0
20	209986	A	T	227	0	0	2	0.008734	229	106	0	0	1
20	210097	C	T	0	82	0	1	0.012048	83	0	37	0	0
...
```

Last four columns represent the number of reads, for each base, that are on the reverse strand. This information can be used to compute strand bias at base-specific resolution.

####  SNPs pileup
Provides pileup information for all positions specified in the input VCF and uses the alternative alleles specified in the VCF file for the VAFs calculations. 

```
chr	pos	rsid	ref	alt	A	C	G	T	af	cov
20	68351	rs757428359	A	G	130	0	0	0	0.000000	130
20	68363	rs200192457	A	T	129	0	0	0	0.000000	129
20	68373	rs745889706	T	C	0	0	0	130	0.000000	130
20	68375	rs754912258	A	G	104	0	0	0	0.000000	104
20	68396	rs138777928	C	T	0	141	0	0	0.000000	141
20	68397	rs748102612	G	A	0	0	141	0	0.000000	141
20	68406	rs771803424	A	G	140	0	0	0	0.000000	140
20	76654	rs564320474	G	T	0	0	31	0	0.000000	31
20	76658	rs745496891	C	A	0	49	0	0	0.000000	49
...
```

## Visual reports
PaCBAM includes a script to generate visual data reports written in python.  
It provides different graphs for every output file:  

	rc: gc content and region coverage distributions  
	snps: total SNPs count, total distribution and quantile distributions of alternative heterozygous and alternative homozygous SNPs 
	snvs: base modification count and strand bias distribution  
	pileup: cumulative coverage and allelic fraction distributions  

#### Requirements
Python 2.7.12  
Numpy 1.14.2  
matplotlib 2.2.2  

#### Usage
The report scripts expects as input the prefix of the output files from PaCBAM and the mode in which it was runned.

```
Usage:
 ./report.py -i/--input string -m/--mode int [-o/--output string] [-s/--strandBias]

-i INPUT, --input INPUT
	Specify the input file prefix
-m MODE, --mode MODE
	Specify the mode used
-o OUTPUT, --output OUTPUT
	Specify the output file name (Default input.pdf)
-s, --strandBias
	Plots the strand bias distribution 
```

Mode option:  
	0 Files: .rc, .snps and .snvs  
	1 Files: .rc, .snps, .snvs and .pileup  
	2 Files: .snps  
	3 Files: .rc  
	4 Files: .pileup  

StrandBias reporting is available only in modes 0 and 1.

#### Example
The following command computes the visual reports for the example data.

```
./report.py -i example/NGSData -m 1 -o reports/reports.pdf

```

#### Output file
The report script produces a single pdf file with all the graphs of the choosen mode.

![cumulativeCoverage](https://bitbucket.org/CibioBCG/pacbam/raw/master/reports/cumulativeCoverage.png)  
*Example of PaCBAM reporting the cumulative coverage distribution for all positions reported in the PaCBAM pileup output file.*

![SNPsTypes](https://bitbucket.org/CibioBCG/pacbam/raw/master/reports/SNPsTypes.png)  
*Example of PaCBAM reporting on allelic fraction (AF) distribution of all positions contained in the PaCBAM SNPs output file. SNPs are classified as heterozygous or alternative homozygous based on standard AF thresholds. Classification is also reported stratified by coverage quartiles.*

![baseModification](https://bitbucket.org/CibioBCG/pacbam/raw/master/reports/baseModification.png)  
*Example of PaCBAM reporting on distribution of alternative bases found for each reference base across all positions reported in the SNVs PaCBAM output file (i.e. all positions with non-zero variant allelic fraction).*

![regionCoverage](https://bitbucket.org/CibioBCG/pacbam/raw/master/reports/regionCoverage.png)  
*Example of PaCBAM reporting on mean depth of coverage distribution computed across all regions reported in the genomic regions of the PaCBAM output file. Distribution is reported both for regions overall mean coverage and for regions fractions maximizing mean coverage.*

## Licence
 
PaCBAM is released under [MIT](https://bitbucket.org/CibioBCG/pacbam/src/master/COPYING) licence.

## Credit
