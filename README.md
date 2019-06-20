

[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/bioinfologics/bsg/blob/master/LICENSE)
[![codecov](https://codecov.io/gh/bioinfologics/bsg/branch/master/graph/badge.svg)](https://codecov.io/gh/bioinfologics/sdg)
[![Build Status](https://travis-ci.org/bioinfologics/sdg.svg?branch=master)](https://travis-ci.org/bioinfologics/sdg)


# SDG
The **S**equence **D**istance **G**raph is a framework to work with genome graphs and raw reads.

It provides a workspace that can contain a graph and datastores for paired, linked and long reads. These reads can be mapped to the graph.

A SWIG API enables SDG to be used as a Python module, and there is experimental Julia and R support.

The documentation is severily lacking right now. We're working on it.

There is [SDG documentation](https://bioinfologics.github.io/sdg/), although again it is lacking right now.

# Using SDG 
## Graph simplification/scaffolding using datastores and existing heuristics

**Warning: SDG is undergoing major changes and this master repository is outdated. I will be updated soon with a working version.**

The current version of needs as input:
- A GFA+fasta graph of the current (input) assembly. Preferably a w2rap-contigger contigs_raw.gfa graph assembled at k=200.
- An illumina dataset to use as k-mer spectrum for the KCI. This is usually the dataset that generated the input assembly graph.
- 10x Chromium Illumina and/or LMP and/or Pacbio and/or Nanopore reads.

All read sequence files need to be converted to datastores. These are
indexed files SDG can use to access every read's sequence
efficiently.

A full workspace must be created, with all of these components, this can later save other
pre-computed status data such as read mapping positions.

Some example commands (**warning, SDG is being thoroughly updated right now, this may only be useful as a trivial exaple now**):

```
#create new workspace with graph
sdg-workspace make -g {input.gfa} -o {wsname}

#create new datastore with linked reads
sdg-datastore make -t 10x -1 {10x_R1.fastq} -2 {10x_R2.fastq} -o {dsname.sdgds}

#add 10x datastore to ws
sdg-workspace add -w {wsname.sdgws} -l {dsname.sdgds}

#add KCI spectrum computed from PE files into workspace
sdg-kmerspectrum make -w {wsname.sdgws} -f {r1.fastq} -f {r2.fastq}

#map reads to graph and save ws to new _mapped file
sdg-mapper -w {wsname.sdgws} -o {wsname_mapped}
```
