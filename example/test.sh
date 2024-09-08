#!/bin/bash

set -e

TESTDATADIR=$(dirname ${0});

# Decompress the test data
pushd $TESTDATADIR;
for FILE in $(ls -1 *gz)
do 
	tar -zxvf ${FILE};
done
popd

mkdir tmp

# Generate the new test data
./pacbam bam=${TESTDATADIR}/NGSData.bam bed=${TESTDATADIR}/TargetRegions.bed vcf=${TESTDATADIR}/SNPsInTargetRegions.vcf fasta=${TESTDATADIR}/human_g1k_v37.fasta mode=1 genotype out=tmp/

# Test the differences
for GZFILE in $(ls -1 tmp/NGSData.*)
do 
	FILE=$(basename $GZFILE .gz);
	diff -q tmp/${FILE} ${TESTDATADIR}/${FILE}
done

# Clean up the testdata
rm ${TESTDATADIR}/human_g1k_v37.fasta*; 

rm -r tmp
