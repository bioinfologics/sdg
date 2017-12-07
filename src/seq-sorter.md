# seq-sorter output files

## Link generation diagram
@image html block_match_diagram.png

### read_block.log

To identify why blocks are created. Each line is a read (in the same order that the input file).


`blockCtgDif` Number of blocks created by jumping to a new contig

`blockOffsetDif` Number of blocks created because a change of offset

`Validblock` Number of valid blocks in the read

`Tblock` Total number of blocks generated on the read

the parameter `--min_kmers_to_call_match uint` is the number of matched kmers required to build a valid block, will affect the number of `Validblock` generated and the output of this file. Increasing `--min_kmers_to_call_match` will decrease the number of blocks but each block will be more robust.


```
blockCtgDif,blockOffsetDif,Validblock,Tblock
2,2,2,3
0,0,1,1
0,0,0,1
0,0,1,1
9,7,1,10
2,2,3,3
0,0,1,1
7,4,3,8
1,1,1,2
12,11,1,13
```

### read_link.log

Links generated per read, this is the number links that a reads generates between different contigs. Each line is a read. For a read to generate a link a read mus contain at least a valid block in two different contigs.

The parameter `--min_kmers_to_call_match` will affect indirectly the number of links through changing the number of valid blocks.

```
Number of links
1
0
0
0
0
2
0
2
0
```

### read_coverage.log

This file show the amount of matched kmers and the kmers that formed a valid block for each read. Each line is a different read.

`Mappedkmers` Number of matched kmers between the read and any contig

`validBlockKmers` Number of matched kmers within all valid blocks

`Tkmers` Total number of kmers in the read (# RL-K+1)

The parameter `--min_kmers_to_call_match` will affect the number of `validBlockKmers ` through changing the number of valid blocks.

```
Mappedkmers,validBlockKmers,Tkmers
321,312,12125
24,24,13050
3,0,1438
698,698,15166
153,142,5071
123,123,9603
15,15,2350
827,790,16549
333,327,11093
```

### read_length.log

Just the length of the reads.

```
Name,Length
SRR6031717.1,12155
SRR6031717.2,13080
SRR6031717.3,1468
SRR6031717.4,15196
SRR6031717.5,5101
SRR6031717.6,9633
SRR6031717.7,2380
SRR6031717.8,16579
SRR6031717.9,11123
```


### read_contig_link.txt

This file contains the links created between contigs. 

A link is formed when there are more than `--min_count` and less than `--max_count` reads that form valid blocks in 2 contigs.

```
Links
-3891;-8197
-4414;6134
-5606;-6134
-2249;10615
10615;-11781
-3454;6946
3130;-6946
2820;3130
57;-797
```

### Final kc

A file containing the links. This file's format is a uint64_t with the number of links on the file. 
Followed by the binary representation of the link objects.

### graph.gfa

Final graph with the connections.