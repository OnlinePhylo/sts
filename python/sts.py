#!/usr/bin/env python
from __future__ import division

from Bio import Phylo
import collections
import itertools
import pydot
import json
import sys

class TooFewPredecessors(ValueError):
    pass

class NoPathError(ValueError):
    pass

def parse_fields_and_data(obj):
    fields = obj['fields']
    return [dict(itertools.izip_longest(fields, item)) for item in obj['data']]

class Edge(object):
    def __init__(self, length, node):
        self.length = length
        self.node = node

    def __repr__(self):
        return '<Edge %#x: %g -> %r>' % (id(self), self.length, self.node)

class Node(object):
    def __init__(self, id, edge1=None, edge2=None, name=None):
        if (edge1 is None) != (edge2 is None):
            raise ValueError('both or neither of edge1 and edge2 must be None')
        self.id = id
        self.edge1 = edge1
        self.edge2 = edge2
        self.name = name

    def __repr__(self):
        return '<Node %s(%#x, %r): %r, %r>' % (
            self.id, id(self), self.name, self.edge1, self.edge2)

class State(object):
    def __init__(self, id, node, predecessor=None):
        self.id = id
        self.node = node
        self.predecessor = predecessor

    def nth_predecessor(self, n):
        ret = self
        for x in xrange(n):
            ret = ret.predecessor
            if ret is None:
                raise TooFewPredecessors()
        return ret

    def __repr__(self):
        return '<State %s(%#x): %r -> %r>' % (
            self.id, id(self), self.predecessor, self.node)

class Generation(object):
    @classmethod
    def of_json_object(cls, state_log, obj):
        self = cls()
        node_map = state_log.node_map
        for node in parse_fields_and_data(obj['nodes']):
            if node['id'] in node_map:
                continue
            edge1 = edge2 = None
            if node['child1'] is not None:
                edge1 = Edge(node['length1'], node_map[node['child1']])
            if node['child2'] is not None:
                edge2 = Edge(node['length2'], node_map[node['child2']])
            node_map[node['id']] = Node(node['id'], edge1, edge2, node['name'])

        state_map = state_log.state_map
        to_resolve = []
        for state in parse_fields_and_data(obj['states']):
            if state['id'] in state_map:
                continue
            state_map[state['id']] = state_obj = State(
                state['id'], node_map.get(state['node']), state['predecessor'])
            to_resolve.append(state_obj)
        for state in to_resolve:
            state.predecessor = state_map.get(state.predecessor)
        self.particles = [state_map[particle] for particle in obj['particles']]
        self.unique_states = set(self.particles)
        self.unique_predecessors = {s.predecessor for s in self.unique_states if s.predecessor is not None}
        return self

    def __repr__(self):
        return '<Generation %#x of %d particles>' % (id(self), len(self.particles))

class StateLog(object):
    @classmethod
    def of_json_file(cls, infile):
        self = cls()
        j = json.load(infile)
        generations = j.pop('generations')
        self.metadata = j
        self.node_map = {}
        self.state_map = {}
        self.generations = [Generation.of_json_object(self, obj) for obj in generations]
        generation_particle_counts = {len(g.particles) for g in self.generations}
        if len(generation_particle_counts) != 1:
            raise ValueError("not all generations have the same number of particles")
        self.n_particles, = generation_particle_counts
        return self

    def __repr__(self):
        return '<StateLog %#x of %s nodes; %d states; %d generations>' % (
            id(self), len(self.node_map), len(self.state_map), len(self.generations))

    def average_survival(self):
        survival_num = survival_denom = 0
        for e, g in enumerate(self.generations[:-1]):
            h = self.generations[e + 1]
            survival_num += sum(1 for s in g.particles if s in h.unique_predecessors)
            survival_denom += len(g.particles)
        return survival_num / survival_denom

    def average_mrca_depth(self):
        last_generation = self.generations[-1]
        ancestor_map = {}
        for base in last_generation.unique_states:
            state, depth = base, 1
            while True:
                ancestor_map.setdefault(depth, {}).setdefault(state, set()).add(base)
                state = state.predecessor
                if state is None:
                    break
                depth += 1

        mrca_num = mrca_denom = 0
        for depth, state_map in sorted(ancestor_map.iteritems()):
            joins = collections.defaultdict(list)
            for state, bases in state_map.iteritems():
                joins[state.predecessor].append(len(bases))
            for join_list in joins.itervalues():
                for a, b in itertools.combinations(join_list, 2):
                    mrca_num += depth * a * b
                    mrca_denom += a * b
        return mrca_num / mrca_denom

class GMLeaf(object):
    is_leaf = True

    def __init__(self, name):
        self.name = name

    def __repr__(self):
        return '<GMLeaf(%r) at %#x>' % (self.name, id(self))

    def __str__(self):
        return self.name

class GMMerge(object):
    is_leaf = True

    def __init__(self, a, b):
        self.merged = a, b

    def __repr__(self):
        return '<GMMerge(%r, %r) at %#x>' % (self.merged + (id(self),))

    def __str__(self):
        return 'merge(%s, %s)' % self.merged

class GMNode(object):
    is_leaf = False

    def __str__(self):
        return 'node(%#x)' % (id(self),)

class GMTree(object):
    def __init__(self):
        self.adjacent_nodes = collections.defaultdict(set)
        self.leaves = set()

    def _add_node_to(self, node, other):
        if node.is_leaf:
            self.leaves.add(node)
        self.adjacent_nodes[other].add(node)

    def add_edge(self, a, b):
        self._add_node_to(a, b)
        self._add_node_to(b, a)

    def _remove_node_from(self, node, other):
        nodes = self.adjacent_nodes[other]
        nodes.discard(node)
        if not nodes:
            if other.is_leaf:
                self.leaves.discard(other)
            del self.adjacent_nodes[other]

    def remove_edge(self, a, b):
        self._remove_node_from(a, b)
        self._remove_node_from(b, a)

    def remove_node(self, node):
        if node.is_leaf:
            self.leaves.discard(node)
        others = self.adjacent_nodes.pop(node)
        for other in others:
            self._remove_node_from(node, other)

    def adjacent_via(self, node, via):
        return self.adjacent_nodes[node] - {via}

    def find_path(self, a, b):
        seen = {a}
        stack = [(a, (a,))]
        while stack:
            cur, path = stack.pop()
            for n in self.adjacent_nodes[cur] - seen:
                seen.add(n)
                new_path = path + (n,)
                if n is b:
                    return new_path
                stack.append((n, new_path))
        raise NoPathError(a, b, 'no path could be found')

    def merge(self, a, b):
        if not (a.is_leaf and b.is_leaf):
            raise ValueError('can only merge leaves')
        to_remove = set(self.find_path(a, b))
        merged = GMMerge(a, b)
        new_adjacent = {merged}
        for node in to_remove:
            new_adjacent.update(self.adjacent_nodes[node])
            self.remove_node(node)
        new_adjacent.difference_update(to_remove)
        new_node = GMNode()
        for node in new_adjacent:
            self.add_edge(node, new_node)
        return merged

    def find_k_distance_merges(self, k):
        if k < 2:
            raise ValueError('k must be >= 2 (was %r)' % (k,))
        some_leaf = next(iter(self.leaves))
        seen = {some_leaf}
        stack = [some_leaf]
        distances = collections.defaultdict(collections.Counter)
        merges = set()
        while stack:
            cur = stack.pop()
            for n in self.adjacent_nodes[cur] - seen:
                for m, d in distances[cur].iteritems():
                    if d >= k:
                        continue
                    distances[n][m] = distances[m][n] = d + 1
                    if n.is_leaf and m.is_leaf and d + 1 == k:
                        merges.add((n, m))
                distances[cur][n] = distances[n][cur] = 1
                seen.add(n)
                stack.append(n)

        return merges

    @classmethod
    def of_newick_file(cls, infile):
        self = cls()
        clade_map = {}
        def get_node(clade):
            node = clade_map.get(clade)
            if node is None:
                if clade.is_terminal():
                    node = GMLeaf(clade.name)
                else:
                    node = GMNode()
                clade_map[clade] = node
            return node

        for clade in Phylo.read(infile, 'newick').find_clades():
            clade_node = get_node(clade)
            for child in clade:
                self.add_edge(clade_node, get_node(child))

        return self

    def to_dot(self):
        graph = pydot.Dot(graph_type='digraph')
        for a, bs in self.adjacent_nodes.iteritems():
            for b in bs:
                graph.add_edge(pydot.Edge(str(a), str(b)))
        return graph

    def clone(self):
        clone = type(self)()
        clone.adjacent_nodes.update(self.adjacent_nodes)
        clone.leaves.update(self.leaves)
        return clone

def main():
    with open(sys.argv[1]) as fobj:
        state_log = StateLog.of_json_file(fobj)

    print 'average survival', state_log.average_survival()
    print 'average mrca depth', state_log.average_mrca_depth()

if __name__ == '__main__':
    main()
