.PHONY: all

all: 4taxon.fasta 4taxon.nex 4taxon.tre mb_4taxon.p mb_4taxon.t mb_4taxon.cherries.csv mb_4taxon.ext.svg

4taxon.tre:
	R --slave -e 'set.seed(1);library(ape);write.tree(rtree(4, br=rexp, rate=5), "4taxon.tre");'

4taxon.fasta: bppseqgen.bpp
	bppseqgen param=$^

4taxon.ML.dnd: bppml.bpp 4taxon.fasta
	bppml param=$^

%.nex: %.fasta
	seqmagick convert --output-format nexus --alphabet dna $< $@

mb_%.t mb_%.p: %.nex 4taxon.mb
	mb ./4taxon.mb

%.cherries.csv: %.t 4taxon.fasta
	../../_build/natural-extension input.sequence.file=4taxon.fasta \
		natural_extension.exp_mean=10.0 \
		natural_extension.burnin=0 \
		natural_extension.trees_nexus=$< \
		natural_extension.output_path=$@

%.ext.svg: plot_cherries.R %.cherries.csv %.p
	Rscript $^ $@
