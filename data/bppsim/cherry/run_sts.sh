GSL_RNG_SEED=1 ../../../_build/sts seq.fasta -p 199 | sort -n -r | uniq > sts.sort.uniq.out
head -n 1 sts.sort.uniq.out | cut -f 2 > sts.tre
