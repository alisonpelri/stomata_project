# run CD-HIT
## experiment different seq identity threshold

for i in `ls *fa | sed 's/\.fa//'`
do
	### 1.0
	cd-hit -i ${i}.fa -o ${i}_100 -c 1 -n 5 -T 5
	### 0.95
	cd-hit -i ${i}.fa -o ${i}_95 -c 0.95 -n 5 -T 5
	#0.90
	cd-hit -i ${i}.fa -o ${i}_90 -c 0.9 -n 5 -T 5
done

# linearize the multifasta files
for i in `ls *fa`; do awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $i > ${i}sta; done
