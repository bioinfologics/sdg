Workspace
=============

The WorkSpace holds the information regarding a project in memory and can be written to and loaded from disk using a file.

It contains:

    - A SequenceDistanceGraph which holds all the sequences, and links between them. This is considered the base graph.
    - Multiple DistanceGraphs they contain alternative links between the sequences in the SDG nodes which can come from one or more sources of information (datastores).
    - Paired, Linked and Long reads datastores which contain all the sequencing information.
    - KmerCounters which contain the kmers used for the count in the form of an index and multiple counts for those kmers. These counts can come from one or more data sources.
    - A Journal which contain details about any operation applied to the WorkSpace that transformed it in any way.

Graph
******

The most important features of this class are:

    - Nodes start from 1.
    - Nodes are signed to represent the direction in which they are being considered.

Short reads
************

Short reads are stored in files which allow fast random access to the sequences.

Linked reads
*************

These are stored in files which are tag sorted to allow random access to both the sequences of the reads and also the sequences of any 10x tag.

Long reads
*************

These are stored in files that hold an index to the location of the records of the reads to allow for fast random access.

Journal
**********

The Journal holds information of any transformative operation that has been applied to a WorkSpace.