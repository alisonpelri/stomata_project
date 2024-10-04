# Stomatal Evolution Pipeline

## Description
This repository contains a series of scripts for analyzing the evolution of stomata in different plant species. The pipeline covers tasks from downloading raw fastq files to transcriptome assembly, translation, and gene tree inference.

### Main Scripts
- `download_fq.sh`: Downloads raw fastq files from SRA and concatenates reads.
- `quality_control.sh`: Performs quality control on the reads using FastQC, MultiQC, and Trimmomatic.
- `transcriptome_assembly.sh`: Assembles the transcriptome using Trinity.
- `translate_orfs.sh`: Translates the longest ORFs from Trinity-assembled transcripts.
- `processing_proteomes.sh`: Removes redundancy from predicted proteomes and linearize sequences 
- `Osativa_homology_search.sh`: Search for homologous sequences on the Orya sativa proteomes using Arabidopsis thaliaana protein sequences as querries
- `homology_search.sh`: Search for homologous sequences in the whole proteomes database using O. sativa sequences as querries
- `gene_tree_inference.sh`: Infers gene trees using mafft, trimal, and iqtree.
- `phylogeny_plots.R`: Each block processes a phylogenetic tree, and you can easily add or modify tree plots for different genes by following the same structure.

### Requirements
- [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [MultiQC](https://multiqc.info/)
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- [TransDecoder](https://github.com/TransDecoder/TransDecoder)
- [mafft](https://mafft.cbrc.jp/alignment/software/)
- [iqtree](http://www.iqtree.org/)
- [CD-HIT](http://weizhongli-lab.org/cd-hit/)
- [BUSCO](https://busco.ezlab.org/)
