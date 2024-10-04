#!/bin/bash
# Author: Alison Menezes
# Description: Infer gene trees using mafft, trimal, and iqtree.

# Gene tree inference and alignment
for i in `ls -Sr *fa | sed 's/_homologs\.fa//'`; do
#align protein sequences
    mafft --localpair --maxiterate 1000 ${i}_homologs.fa > ${i}.aln
#trim alignments
    trimal -in ${i}.aln -out ${i}_trimal.fa -gappyout
#phylogeny reconstruction
    iqtree -s ${i}_trimal.fa -B 1000 -m TESTMERGE -T AUTO --threads-max 20 --prefix $i --redo
done
