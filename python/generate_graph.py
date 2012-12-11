#!/usr/bin/env python
from __future__ import division

import argparse
import collections
import sys

from analyze import StateLog
import pydot

def main():
    p = argparse.ArgumentParser()
    p.add_argument('json_log', type=argparse.FileType('r'))
    p.add_argument('output_png')
    a = p.parse_args()

    with a.json_log  as fobj:
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

    graph.write_png(a.output_png)

main()
