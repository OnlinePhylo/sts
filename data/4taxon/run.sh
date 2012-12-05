#!/bin/sh

set -e
set -u

R --slave -e 'set.seed(1);library(ape);write.tree(rtree(4, br=rexp, rate=5), "4taxon.tre");'
bppseqgen param=bppseqgen.bpp
bppml param=bppml.bpp
seqmagick convert --output-format nexus --alphabet dna 4taxon.fasta 4taxon.nex
mb ./4taxon.mb
