[![codecov](https://codecov.io/gh/ljyanesm/sg/branch/master/graph/badge.svg)](https://codecov.io/gh/ljyanesm/sg)
[![Build Status](https://travis-ci.org/ljyanesm/sg.svg?branch=master)](https://travis-ci.org/ljyanesm/sg)
# bsg
A basic sequence graph, with some associated tools.

Information on the tools on this file, information on sglib can be found
on the [sglib readme](src/sglib/README.md)


# Running the untangler (graph simplification)

The current version of needs as input:
- A GFA+fasta graph of the current (input) assembly. Preferably a
w2rap-contigger contigs_raw.gfa graph assembled at k=200.
- An illumina dataset to use as k-mer spectrum for the KCI. This is
usually the datasdt that generated the input assembly graph.
- 10x Chromium Illumina.

All read sequence files need to be converted to datastores. These are
indexed files BSG can use to access every read's sequence
efficiently.

A full bsg-scaffolder workspace must be created, with all of these
components, this can later save other
pre-computed status data such as read mapping positions.

Commands to run untangler2 with default settings for 10x
(WARNING: the commands do not work like this yet!)

```
#create new workspace with graph
bsg-workspace make -g {input.gfa} -o {wsname}

#create new datastore with linked reads
bsg-datastore make -t 10x -1 {10x_R1.fastq} -2 {10x_R2.fastq} -o {dsname.bsgds}

#add 10x datastore to ws
bsg-workspace add -w {wsname.bsgws} -l {dsname.bsgds}

#add KCI spectrum computed from PE files into workspace
bsg-kmerspectrum make -w {wsname.bsgws} -f {r1.fastq} -f {r2.fastq}

#map reads to graph and save ws to new _mapped file
bsg-mapper -w {wsname.bsgws} -o {wsname_mapped}

#create flows and save ws to new _flows file
bsg-flowmaker -w {wsname_mapped.bsgws} -o {wsname_flows}

#finally, run untangler
bsg-untangler -w {wsname_flows.bsgws} -o {untangled}
```


%## Format requirements for the GFA+fasta input graph

%## Generating a read kmer spectrum for KCI

%## Generating datastores from fastq files

%## Creating a workspace

%## Launching a scaffolding run