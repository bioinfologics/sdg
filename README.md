

[![](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/bioinfologics/bsg/blob/master/LICENSE)
[![](https://codecov.io/gh/bioinfologics/sdg/branch/master/graph/badge.svg)](https://codecov.io/gh/bioinfologics/sdg)
[![](https://travis-ci.org/bioinfologics/sdg.svg?branch=master)](https://travis-ci.org/bioinfologics/sdg)


# SDG
The **S**equence **D**istance **G**raph is a framework to work with genome graphs and sequencing data. It provides a workspace build around a Sequence Distance Graph, datastores for paired, linked and long reads, read mappers, and k-mer counters. It can be used to perform different types of sequence analyses.

SDG can be used as a Python module through its SWIG API. R and Julia support are experimental and should not be used at this stage.

The documentation is severily lacking right now. We're working on it. There is [SDG API documentation](https://bioinfologics.github.io/sdg/), although again it is lacking right now.

# WARNING

**SDG is undergoing major changes right now preparing for v1. Please come back in a few weeks for serious usage.**

# Installation



# Usage



## Command line tools

### sdg-datastore

Tool to create datastores from raw reads. It also creates a special type of datastore: the KmerCounts, which computes frequency of k-mers from a read set. Since the k-mers are only counted if they appear in a WorkSpace's SequenceDistanceGraph, you will need a workspace to create a KmerCounts.

### sdg-workspace

Tool to create and modify WorSpaces. The WorkSpaces bind together a SequenceDistanceGraph, aditional DistanceGraphs defined over its nodes, read datastores (and their mappings) and k-mer count datastores. You will normally use WorkSpaces to save intermediate status of pipelines and such.

### sdg-mapper

A simple tool to map reads to the SequenceDistanceGraph within a Workspace.

### sdg-dbg

This tool takes a read datastore as input and constructs a DBG, then maps the reads back to it, creates a k-mer count and writes down a neat WorkSpace with all this.



## Python module examples

### Example #1: short and long read hybrid assembly/scaffolding

For this example, you need to first have both the illumina and pacbio reads in datastores, and then create a starting WorkSpace with the DBG of the read sequences:

```bash
sdg-datastore make -t paired -o ecoli_pe -1 ../ecoli_pe_r1.fastq -2 ../ecoli_pe_r2.fastq -s 301
sdg-datastore make -t long -o ecoli_pb -L ../ecoli_pb_all.fastq
sdg-dbg -p ecoli_pe.prseq -o ecoli_assm
```

Now the long reads datastore can be added, the reads mapped, the linkage created, and the repeats eliminated, to create a scaffolded assembly that will contain only the non-repetitive nodes with mappings from the original DBG:

```python
import pysdg as sdg

# Load sdg-dbg's workspace from disk, add the pacbio datastore
ws=sdg.WorkSpace('ecoli_assm.bsgws')
ws.long_read_datastores.push_back(sdg.LongReadsDatastore(ws,'ecoli_pb.loseq'))

# Map long reads, filter by 10% id, pick best matches and rechain
lds=ws.long_read_datastores[0]
lds.mapper.k=11
lds.mapper.map_reads()
lds.mapper.filter_mappings_by_size_and_id(500,10)
lds.mapper.improve_filtered_mappings()

# Create a LinkageUntangler, with linkage from the mapped long reads
u=sdg.LinkageUntangler(ws)
lr_mldg=u.make_longRead_multilinkage(lds.mapper)

# Select large nodes, any CI, create linkage between 1st neighbours on selection
u.select_nodes_by_size_and_ci(1100,0,100)
nsl=u.make_nextselected_linkage(lr_mldg)
nsl.write_to_gfa('lr_scaffolded_with_repeats.gfa')

# Deselect nodes with many inputs or outputs (repeats), create linkage with no repeats
for n in range(len(ws.sdg.nodes)):
    if u.selected_nodes[n] and ( len(nsl.get_bw_links(n))>1 or len(nsl.get_fw_links(n))>1):
      u.selected_nodes[n]=False
nsl_nr=u.make_nextselected_linkage(lr_mldg)
nsl_nr.write_to_gfa1('lr_scaffolded_no_repeats.gfa')
```


### Example #2: phasing a trio child genome using k-mer counts

