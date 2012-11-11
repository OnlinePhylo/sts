GSL_RNG_SEED=1 ../../../_build/sts combined.fasta -o sts.out
sort -n -r sts.out | uniq > sts.sort.uniq.out
head -n 1 sts.sort.uniq.out | cut -f 2 > sts_on_combo.tre
