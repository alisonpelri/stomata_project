#!/bin/bash
# Author: Alison Menezes
# Description: Download fastq files from SRA and concatenate them.

# Downloading reads from SRA
# -p shows progress
for i in `ls -d *`; do
    for j in `cat ${i}/${i}.ids`; do
        printf "dumping $j to $i\n"
        fasterq-dump ${j} -O ${i} -e 20 -p
    done
done

# Concatenating SRP192700 and SRP388644 
## SRP192700
cat SRR89064*1.fastq > SRP192700_R1.fq
cat SRR89064*2.fastq > SRP192700_R2.fq
## SRP388644
cat SRR2070061* > SRP388644.fq
## SRP105098
for i in `ls *_1.fastq | sed 's/_1\.fastq//'`; do
    cat ${i}_1.fastq >> SRP105098_R1.fq
    rm ${i}_1.fastq
    cat ${i}_2.fastq >> SRP105098_R2.fq
    rm ${i}_2.fastq
done
