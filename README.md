# SharpSAT-TD

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4880703.svg)](https://doi.org/10.5281/zenodo.4880703)

Submission to model counting competition 2021 by Tuukka Korhonen and Matti Järvisalo (University of Helsinki).
SharpSAT-TD is based on [SharpSAT](https://github.com/marcthurley/sharpSAT), with the main new features being the use of tree decompositions in decision heuristics, new preprocessor, and directly supporting weighted model counting.


SharpSAT-TD supports exact model counting, exact weighted model counting with arbitrary precision floats, and exact weighted model counting with doubles.
See a detailed description in [description.pdf](https://github.com/Laakeri/sharpsat-td/blob/main/description.pdf).


# Compiling

The external dependencies needed are the [GMP library](https://gmplib.org/), the [MPFR library](https://www.mpfr.org/), and [CMAKE](https://cmake.org/).

To compile and link dynamically use

``./setupdev.sh``

To compile and link statically use

``./setupdev.sh static``


The binaries sharpSAT and flow_cutter_pace17 will be copied to the [bin/](https://github.com/Laakeri/sharpsat-td/tree/main/bin) directory.

# Running

The currently supported input/output formats are those of [Model counting competition 2021](https://mccompetition.org/assets/files/2021/competition2021.pdf).


Example unweighted model counting:
`cd bin`
`./sharpSAT -decot 1 -decow 100 -tmpdir . -cs 3500 ../examples/track1_009.cnf`


Example weighted model counting with arbitrary precision:
`cd bin`
`./sharpSAT -WE -decot 1 -decow 100 -tmpdir . -cs 3500 -prec 20 ../examples/track2_003.wcnf`


Example weighted model counting with double precision:
`cd bin`
`./sharpSAT -WD -decot 1 -decow 100 -tmpdir . -cs 3500 ../examples/track2_003.wcnf`


In the competition setting the value of the `-decot` flag was 120.

## Flags

- `-decot` - the number of seconds to run flowcutter to find a tree decomposition. Required. Recommended value 60-600 if running with a total time budjet of 1800-3600 seconds.
- `-tpmdir` - the directory to store temporary files for running flowcutter. Required.
- `-decow` - the weight of the tree decomposition in the decision heuristic. Recommended value >1 if the heuristic should care about the tree decomposition.
- `-cs` - limit of the cache size. If the memory upper bound is X megabytes, then the value here should be around x/2-500.
- `-WE` - enable weighted model counting with arbitrary precision.
- `-WD` - enable weighted model counting with double precision.
- `-prec` - the number of digits in output of weighted model counting. Does not affect the internal precision.