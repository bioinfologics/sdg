.. image:: https://img.shields.io/badge/license-MIT-green.svg
    :target: https://github.com/bioinfologics/bsg/blob/master/LICENSE
.. image:: https://codecov.io/gh/bioinfologics/sdg/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/bioinfologics/sdg
.. image:: https://travis-ci.org/bioinfologics/sdg.svg?branch=master
    :target: https://travis-ci.org/bioinfologics/sdg

Sequence Distance Graph
========================

The Sequence Distance Graph (**SDG**) is a framework to work with genome graphs and sequencing data. It provides a workspace built around a Sequence Distance Graph, datastores for paired, linked and long reads, read mappers, and k-mer counters. It can be used to perform different types of sequence analyses.

SDG can be used as a Python module through its SWIG API. R and Julia support are experimental.

.. warning:: SDG is undergoing major changes right now preparing for v1. Please come back in a few weeks.


Installation
#############

The installation process requires the following dependencies:

- GCC >= 6
- CMake >= 3.14.3
- SWIG == 4.0.0 (To compile the Python interface)
- Python >= 3.7 (To compile the Python interface)
- Git LFS (To obtain the test datasets)

The installation process consists of generating the configuration files using CMake and using make.

.. code-block:: bash

    git pull https://github.com/bioinfologics/sdg
    cd sdg
    mkdir build
    cd build
    cmake -DBUILD_PYTHON_INTERFACE=ON ../
    make
    make test
    make install


Usage
#####

Command line tools
########################


sdg-datastore
*************************

Tool to create datastores from raw reads. It also creates a special type of datastore: the KmerCounts, which computes frequency of k-mers from a read set. Since the k-mers are only counted if they appear in a WorkSpace's SequenceDistanceGraph, you will need a workspace to create a KmerCounts.

sdg-workspace
*************************

Tool to create and modify WorkSpaces. The WorkSpaces bind together a SequenceDistanceGraph, additional DistanceGraphs defined over its nodes, read datastores (and their mappings) and k-mer count datastores. You will normally use WorkSpaces to save intermediate status of pipelines and such.

sdg-mapper
*************************

A simple tool to map reads to the SequenceDistanceGraph within a Workspace.

sdg-dbg
*************************

This tool takes a read datastore as input and constructs a DBG, then maps the reads back to it, creates a k-mer count and writes down a neat WorkSpace with all this.



Python module examples
###########################

Example #1: short and long read hybrid assembly/scaffolding
****************************************************************

For this example, you need to first have both the Illumina and Pacbio reads in datastores, and then create a starting WorkSpace with the DBG of the read sequences:

.. code-block:: bash

    sdg-datastore make -t paired -o ecoli_pe -1 ../ecoli_pe_r1.fastq -2 ../ecoli_pe_r2.fastq -s 301
    sdg-datastore make -t long -o ecoli_pb -L ../ecoli_pb_all.fastq
    sdg-dbg -p ecoli_pe.prseq -o ecoli_assm

Now the long reads datastore can be added, the reads mapped, the linkage created, and the repeats eliminated, to create a scaffolded assembly that will contain only the non-repetitive nodes with mappings from the original DBG:

.. code-block:: python
    :linenos:

    import pysdg as sdg

    # Load sdg-dbg's workspace from disk, add the pacbio datastore
    ws=sdg.WorkSpace('ecoli_assm.sdgws')
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

Example #2: phasing a trio child genome using k-mer counts
*****************************************************************
