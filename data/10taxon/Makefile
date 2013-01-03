.PHONY: all

NATURAL_EXTENSION ?= ../../_build/natural-extension

all: 10taxon.fasta 10taxon.nex 10taxon.tre mb_10taxon.p mb_10taxon.t mb_10taxon.nat_ext.csv mb_10taxon.ext.svg

10taxon.tre:
	R --slave -e 'set.seed(1);library(ape);write.tree(rtree(10, br=rexp, rate=5), "10taxon.tre");'

10taxon.fasta: bppseqgen.bpp 10taxon.tre
	bppseqgen param=$^

10taxon.ML.dnd: bppml.bpp 10taxon.fasta
	bppml param=$^

%.nex: %.fasta
	seqmagick convert --output-format nexus --alphabet dna $< $@

mb_%.p mb_%.t: %.nex 10taxon.mb
	mb ./10taxon.mb

%.nat_ext.csv: %.t 10taxon.fasta $(NATURAL_EXTENSION)
	$(NATURAL_EXTENSION) \
		input.sequence.file=10taxon.fasta \
		natural_extension.exp_mean=10.0 \
		natural_extension.seed=1 \
		natural_extension.burnin=1250 \
		natural_extension.trees_nexus=$< \
		natural_extension.output_path=$@

%.ext.svg: plot_natural_extension.R %.nat_ext.csv %.p
	Rscript $^ $@
