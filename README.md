Sequential Tree Sampler (STS)
=============================

## Dependencies

* [smctc][smctc] - included as git submodule (`git submodule update --init`)
* [beagle][beagle] version 2.1
* [Bio++ version 2.2.0][bpp] `core`, `seq`, and `phyl` modules (Debian packages `libbpp-core-dev libbpp-seq-dev libbpp-phyl-dev`)
* [cmake][cmake]

## Compiling

1. Install dependencies
1. run `make`

Binaries will be build in `_build/release`

## Adding taxa to an existing posterior

The tool `sts-online` adds taxa to an existing posterior tree sample.
`sts-online` operates on a fasta file and tree file in nexus format.
The fasta file must contain an alignment with a superset of the taxa in the tree file.

Currently only the Jukes-Cantor model is supported.

### Example invocation

    _build/release/sts-online 50taxon-01.fasta 50tax_trim.run1.t 50tax_trim.sts.json

Notes
-----
[Talk PDF from Liangliang Wang on Combinatorial SMC][csmc]



[smctc]: http://www2.warwick.ac.uk/fac/sci/statistics/staff/academic-research/johansen/smctc/
[beagle]: https://code.google.com/p/beagle-lib/
[bpp]: http://biopp.univ-montp2.fr/
[csmc]: http://www2.warwick.ac.uk/fac/sci/statistics/crism/workshops/sequentialmontecarlo/programme/smc2012_lwpdf.pdf
[cmake]: http://www.cmake.org/
