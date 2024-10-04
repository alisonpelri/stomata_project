# Translating longest ORFs assembled by Trinity
## extract longest ORFs from Trinity output
for i in `ls -Sr */*_trinity.fasta | sed 's/_trinity\.fasta//'`; do perl get_longest_isoform_seq_per_trinity_gene.pl ${i}_trinity.fasta > ${i}_trinity.pep; done

## identify all long ORFs ( > 80 nt)
for i in `ls -Sr */*_trinity.fasta | sed 's/_trinity\.fasta//'`; do TransDecoder.LongOrfs -t ${i}_trinity.pep -m 80 -O ${i}_transdecoder/; done

## make db with the stomatal development genes from A. thalina and O. sativa
makeblastdb -in Athaliana_447_Araport11.protein_primaryTranscriptOnly.fa -dbtype prot

## run blastp to enable homology-based coding region identification
for i in `ls -Sr */*_trinity.fasta | sed 's/_trinity\.fasta//'`; do blastp -query ${i}_translated.fa -db ath.fa -evalue 1e-5 -outfmt 6 -max_target_seqs 1 -num_threads 10 -out ${i}_blastp.tsv; done

## predict which ORFs are likely to be coding
for i in `ls -Sr */*_trinity.fasta | sed 's/_trinity\.fasta//'`; do TransDecoder.Predict -t ${i}_trinity.fasta -O ${i}_transdecoder/ --retain_blastp_hits ${i}_blastp.tsv --single_best_only; done
