#!/bin/bash
# Author: Alison Menezes
# Description: Perform quality control using FastQC and Trimmomatic.

# Run FastQC
for i in `ls -d SRP*`; do
    fastqc -o ${i}/ ${i}/*fastq
done

# Run multiqc to aggregate FastQC results
multiqc qual_check/fastqc_dir

# Run Trimmomatic to trim paired-end reads
for i in `ls *_1.fastq | sed 's/_1\.fastq//'`; do
    trimmomatic PE -threads 15 ${i}_1.fastq ${i}_2.fastq ${i}_trimmed_R1.fq ${i}_unpaired_R1.fq.gz ${i}_trimmed_R2.fq ${i}_unpaired_R2.fq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50
done
