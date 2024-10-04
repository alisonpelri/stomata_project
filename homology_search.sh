## make db with the stomatal development genes from A. thalina and O. sativa
cat *fa > proteomes_db.fa
makeblastdb -in proteomes_db.fa -dbtype prot

## run blastp search on translated proteomes db
for i in `ls *.fa | sed s/\.fa// | sed 's/Osa//'`; do blastp -query Osa${i}.fa -db translated_proteomes_db.fa -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen" -num_threads 20 -out ${i}_blastp.tsv; done

# Parse BLASTp output
# getting target sequence ids from blastp search only hits with:
## e-value < 1e-10
for i in `ls *_blastp.tsv | sed 's/_blastp\.tsv//'`; do awk '$11 <= 1e-10 {print}' ${i}_blastp.tsv | cut -f 2 | sort | uniq > ${i}.ids; done

# recover protein sequences from predicted proteomes
## use sequences ids to recover sequences from predicted proteomes db
for i in `ls *.ids | sed -e 's/\.ids//'`; do grep -A1 -f ${i}.ids proteomes_db.fa | sed -e '/\-\-/d' > .${i}_homologs.fa; done
