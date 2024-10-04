# make database
## /home/cozzolino/alison/stomata/data/proteomes/predicted/Osativa_db/Osa_db.fa
makeblastdb -in Osa_db.fa -dbtype prot

# blastp search
for i in `ls *.fa | sed -e s'/\.fa//' -e 's/At//'`; do blastp -query At${i}.fa -db Osa_db.fa -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen" -num_threads 5 -out Osa${i}_blastp.tsv; done

# getting target sequence ids from blastp search only hits with:
## e-value < 0.01
## alignment length and percentage of identical positions > 50%
for i in `ls *_blastp.tsv | sed 's/_blastp\.tsv//'`; do awk '$11 <= 1e-10 {print}' ${i}_blastp.tsv | cut -f 2 | sort | uniq > ${i}.ids; done

# recover protein sequences from predicted proteomes
## use sequences ids to recover sequences from Oryza proteome
for i in `ls *.ids | sed -e 's/\.ids//'`; do grep -A1 -f ${i}.ids Osa_db.fa | sed -e '/\-\-/d' > ${i}_homologs.fa; done

# gene tree inference
for i in `ls -Sr *fa | sed 's/_homologs\.fa//'`
do
	mafft --localpair --maxiterate 1000 ${i}_homologs.fa > ${i}.aln
#trim alignments
	trimal -in ${i}.aln -out ${i}_trimal.fa -gappyout
#run iqtree
	iqtree -s ${i}_trimal.fa -B 1000 -m TESTMERGE -T AUTO --threads-max 20 --prefix $i --redo
done
