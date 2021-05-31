# SharpSAT-TD

Submission to model counting competition 2021, unweighted and weighted tracks by Tuukka Korhonen and Matti Järvisalo (University of Helsinki).
SharpSAT-TD is based on [SharpSAT](https://github.com/marcthurley/sharpSAT), with the main new features being the use of tree decompositions in decision heuristics, new preprocessor, and directly supporting weighted model counting.
See a detailed description in [description.pdf](https://github.com/Laakeri/sharpsat-td/blob/main/description.pdf).

# Compiling

The only external dependency needed is the [GMP library](https://gmplib.org/), which can be installed by ``sudo apt-get install libgmp3-dev``.

To compile and link dynamically use
./setupdev.sh
To compile and link statically use
./setupdev.sh static

Everything needed for running will be copied to the [bin/](https://github.com/Laakeri/sharpsat-td/tree/main/bin) directory.

In particular, for unweighted counting, the following are needed:
bin/unweighted.sh
bin/sharpSATU
bin/flow_cutter_pace17

For weighted counting, the following are needed:
bin/weighted.sh
bin/sharpSATW
bin/flow_cutter_pace17

# Running

The main scripts are the bin/unweighted.sh and bin/weighted.sh.

The parameters of these scripts are as documented in the Section 4.4 of [competition2021.pdf](https://mccompetition.org/assets/files/2021/competition2021.pdf).
The parameters --tmpdir and --maxrss are required, others are ignored.