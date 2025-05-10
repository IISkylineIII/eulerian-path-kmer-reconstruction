# Eulerian Path-Based K-mer String Reconstruction

This repository contains a Python script that reconstructs a linear string from a set of paired k-mers by finding an Eulerian path in a directed graph. This method is widely used in bioinformatics, particularly in genome assembly tasks.

## Problem Overview

The problem involves reconstructing a sequence from pairs of overlapping k-mers. These paired k-mers are treated as edges in a directed graph, where each k-mer pair connects a prefix (the first k-mer) to a suffix (the second k-mer). The task is to find an Eulerian path through the graph, which can then be used to reconstruct the original sequence.

## Approach

1. **Graph Construction**: Paired k-mers are used to create a directed graph where each k-mer is a node, and the edges represent the transitions from the prefix to the suffix.
2. **Eulerian Path**: The Eulerian path is found using a depth-first search (DFS) approach, which traverses the graph from a node with an outdegree greater than indegree to reconstruct the original string.
3. **Reconstruction**: The sequence is reconstructed by merging the overlapping k-mers.

## Script

```python
from collections import defaultdict

def string_reconstruction_from_read_pairs(paired_kmers):
    def construct_graph(paired_kmers):
        graph = defaultdict(list)
        for pair in paired_kmers:
            prefix, suffix = pair[0], pair[1]
            graph[prefix].append(suffix)
        return graph

    def find_eulerian_path(graph, start_node):
        path = []
        stack = [start_node]

        while stack:
            current_node = stack[-1]

            if graph[current_node]:
                next_node = graph[current_node].pop()
                stack.append(next_node)
            else:
                path.append(stack.pop())

        return path[::-1]

    graph = construct_graph(paired_kmers)

    # Find a node with an indegree less than outdegree to start the path
    start_node = None
    for node in graph:
        indegree = sum(len(graph[parent]) for parent in graph)
        outdegree = len(graph[node])
        if indegree < outdegree:
            start_node = node
            break

    # If no such node is found, just pick any node as the start
    if start_node is None:
        start_node = list(graph.keys())[0]

    eulerian_path = find_eulerian_path(graph, start_node)

    # Reconstruct the string from the Eulerian path
    reconstructed_string = eulerian_path[0]
    for node in eulerian_path[1:]:
        reconstructed_string += node[-1]

    return reconstructed_string

# Example usage:
paired_kmers = [
    ('ACC', 'ATA'),
    ('ACT', 'ATT'),
    ('ATA', 'TGA'),
    ('ATT', 'TGA'),
    ('CAC', 'GAT'),
    ('CCG', 'TAC'),
    ('CGA', 'ACT'),
    ('CTG', 'AGC'),
    ('CTG', 'TTC'),
    ('GAA', 'CTT'),
    ('GAT', 'CTG'),
    ('GAT', 'CTG'),
    ('TAC', 'GAT'),
    ('TCT', 'AAG'),
    ('TGA', 'GCT'),
    ('TGA', 'TCT'),
    ('TTC', 'GAA')
]

result = string_reconstruction_from_read_pairs(paired_kmers)
print(result)

Example Output
ACCAATTG
Requirements
Python 3.6 or higher

No external libraries are required.

Applications
Genome assembly and sequence reconstruction

Bioinformatics tasks involving paired k-mers

Eulerian path applications in graph theory

License
This project is released under the MIT License.




