#!/usr/bin/env python
from __future__ import division

from analyze import StateLog
import collections
import pydot
import sys

def main():
    with open(sys.argv[1]) as fobj:
        state_log = StateLog.of_json_file(fobj)

    graph = pydot.Dot(graph_type='digraph', rankdir='LR', outputorder='edgesfirst')
    nodes = {}
    for e, g in enumerate(state_log.generations):
        counts = collections.Counter(g.particles)
        for s in g.unique_states:
            nodes[s] = node = pydot.Node('%s x%d' % (s.id, counts[s]), rank=e, shape='point')
            graph.add_node(node)

    for e, g in enumerate(state_log.generations[:0:-1]):
        subgraph = pydot.Subgraph(e, bgcolor='#ff0000')
        graph.add_subgraph(subgraph)
        for s in g.unique_states:
            subgraph.add_edge(pydot.Edge(nodes[s.predecessor], nodes[s]))

    graph.write_png(sys.argv[2])

main()
