# Long read mapper tool

Maps long reads using a minSketch index of the graph generating a list 
of mappings associated with each long read.

## Mapping process

### Generate minSketch index

For each node in the graph generate all minSketches and store them in a map
based on the hash of the minSketch along with it's node and position on the graph.

The resulting index will contain for every minSketch all the node,position pairs in the graph.

### Mapping a single read

All the minSketch of the read are generated and looked for in the graph minSketch index,
if the minSketch matches the index, a mapping is augmented if the ±node is the same as 
the previous minSketch a new mapping is created if the ±node doesn't match.

After all matches are found, the matches are filtered using overlapping window votings.
The window winner extends a block if the previous winner was the same ±node, otherwise
it creates a new block.

Finally, all blocks for a read are stored as ReadMappings.