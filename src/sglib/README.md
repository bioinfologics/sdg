# sglib
A basic sequence graph library.

## Design choices

* Nodes are either a SequenceNode or a SubGraphNode, and directional
* Nodes have a soft-delete is_current flag, there are purge routines
* Sequence nodes by default are canonical, can be expanded into ghost-RC
* Nodes hold information about operations and relationships with other nodes
* Links are directional too, and can define a distance.
* Property classes must implement interface to have IDs updated and be added
to a subscribed list of properties, they will get notified of id updates

## Methods highlight

* read_gfa / write_gfa
* construct_from_unitigs
* add_node / add_ink
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

##

##

## Using the graph

### walking the graph
