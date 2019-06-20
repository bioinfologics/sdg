# sdglib
A sequence distance graph library.

## Design choices

### Nodes and links

We have opted for a representation on which nodes are sequences and go from
sink (-) to source (+) in a canonical representation.

* Nodes are directional, have a left/start/sink side marked '-' and an
end/source side marked '+'.
* Nodes are canonical or palindromic sequences, else they are reverted.
* Links are basically a pair of nodes with signs. +node denotes a connection using
a node's source (end) and -node denotes a connection using a node's sink (start)
* Both nodes and links can be active or historical, and can be purged. See Operations.

### Graph operations

All of this operations will be saved in a journal/log to provide a full
backtrack of the assembly.

* Add node
* Remove node
* Add link
* Remove link
* Edit node: changes the consensus for a node's sequence
* Unfold node: creates extra copies of a node, distributes its links
* Fold nodes: collapses multiple copies of a node, combining their links.
* Split node: breaks a node, creating many smaller parts.
* Combine nodes: merges multiple nodes, creating a single, larger, node

## Methods highlight

* read_gfa / write_gfa
* construct_from_unitigs
* add_node / add_link
* remove_node / remove_link
* purge_node / purge_link / purge_all_nodes / purge_all_links
* sequence_from_path
* single_direction_components: finds single direction components, if a
node is specified the
* replace_path: can be used to expand haplotypes, checks in-out
* solve_region: receives a list of paths through a region, replaces
the region.


## SequenceNode

`SequenceNode n` variables:

* `n.seq`: the sequence, as a string, canonical.
* `n.status`: as of now, always 1, meaning it is active.

## Using the graph

### walking the graph
