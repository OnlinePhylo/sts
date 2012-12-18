#!/usr/bin/env python
from __future__ import division

import collections
import itertools
import json
import sys

class TooFewPredecessors(ValueError):
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

def main():
    with open(sys.argv[1]) as fobj:
        state_log = StateLog.of_json_file(fobj)

    print 'average survival', state_log.average_survival()
    print 'average mrca depth', state_log.average_mrca_depth()

if __name__ == '__main__':
    main()
