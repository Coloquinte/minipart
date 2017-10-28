# Minipart

Minipart is a solver for the hypergraph partitioning problem.

Hypergraph partitioning is a well-known NP-hard optimization problem.
It is an important problem for many fields in computer engineering, like electronic circuit design and task allocation in datacenters.
Minipart is an attempt to bring new algorithms to this field.

Minipart is still in its infancy: it is under active development, focused on solution quality and runtime performance.

## Usage

You can compile Minipart with a recent C++ compiler, CMake and the Boost C++ libraries.

To run it:

    minipart -g <FILE> -o <OUTPUT>

It accepts the [hMetis hypergraph format](http://glaros.dtc.umn.edu/gkhome/fetch/sw/hmetis/manual.pdf).

For more information about the options:

    minipart --help

For example, to run one of the [ISPD98 benchmarks](http://vlsicad.ucsd.edu/UCLAWeb/cheese/ispd98.html) at [2% margin](http://vlsicad.ucsd.edu/GSRC/bookshelf/Slots/Partitioning/HMetisML02Tab.html):

    minipart -g ibm01.hgr --margin 2 

## Name

Minipart is named after [Minisat](http://minisat.se/), the well-known SAT solver. I am not related in any ways to Minisat's authors, but I learned a lot through their codebase and research.

For another open-source partitioner, see [MLPart](http://vlsicad.eecs.umich.edu/BK/PDtools/)
