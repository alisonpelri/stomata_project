#!/bin/bash
# Author: Alison Menezes
# Description: Perform transcriptome assembly using Trinity.

# Run Trinity for paired-end reads
for i in `ls -Sr *R1.fq | sed 's/_R1\.fq//'`; do
    Trinity --seqType fq --left ${i}_R1.fq --right ${i}_R2.fq --CPU 20 --max_memory 100G --full_cleanup --output ${i}_trinity
    if [ -f ${i}_trinity.Trinity.fasta ]; then
        rm ${i}_R1.fq ${i}_R2.fq
    fi
done
