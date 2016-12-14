Sequential Tree Sampler (STS)
=============================

## Dependencies

* [smctc][smctc] - included as git submodule (`git submodule update --init`)
* [lcfit][lcfit] - included as git submodule (`git submodule update --init`)
* [beagle][beagle] version 2.1
* [Bio++ version 2.2.0][bpp] `core`, `seq`, and `phyl` modules. Note that debian & ubuntu up to 16.04 include v2.1.0 which is too old. Bio++ should be installed from source using the `bpp-setup.sh` script on these systems.
* [cmake][cmake]
* [gsl version 1.16][gsl] Note that gsl v2 is not currently supported.
* [nlopt][nlopt]
* [gtest][google test] this is libgtest on debian/ubuntu

## Compiling

1. Install dependencies
2. run `make`

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
[lcfit]: http://github.com/matsengrp/lcfit/
[beagle]: https://code.google.com/p/beagle-lib/
[bpp]: http://biopp.univ-montp2.fr/
[csmc]: http://www2.warwick.ac.uk/fac/sci/statistics/crism/workshops/sequentialmontecarlo/programme/smc2012_lwpdf.pdf
[cmake]: http://www.cmake.org/
[gsl]: https://www.gnu.org/software/gsl/
[nlopt]: http://ab-initio.mit.edu/wiki/index.php/NLopt
[gtest]: https://github.com/google/googletest

