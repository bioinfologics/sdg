# bsg-scaffolder readme


## Requirements for running the scaffolder

The minimum requirements for running the scaffolder are:
- An input graph in GFA+fasta format.
- An input kmer spectrum to be used by the KCI.
- Datastores with reads for scaffolding (at least one of):
    - 10x Chromium Illumina PE reads.
    - Illumina paired reads (already pre-processed if LMP).
    - Pacbio/Nanopore Long reads.

All read sequence files need to be converted to datastores. These are
indexed files the scaffolder can use to access every read's sequence
efficiently.

A full bsg-scaffolder workspace must be created, with all of these
components, this can later save other
pre-computed status data such as read mapping positions.

## Format requirements for the GFA+fasta input graph

## Generating a read kmer spectrum for KCI

## Generating datastores from fastq files

## Creating a workspace


## Launching a scaffolding run

Every scaffolding run will:

1) Load a workspace.
2) Execute the heuristics as asked.
3) Dump a final workspace or graph.
