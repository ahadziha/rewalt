"""
Implements string diagram visualisations.
"""

import networkx as nx

from rewal import utils
from rewal.ogposets import (GrSet, GrSubset)
from rewal.diagrams import Diagram


class StringDiagram:
    """
    Class for string diagrams.
    """
    def __init__(self, diagram, **kwargs):
        utils.typecheck(diagram, {'type': Diagram})

        dim = diagram.dim
        base = nx.DiGraph()
        if dim >= 0:
            for x in diagram.shape[dim]:
                base.add_node(
                        x,
                        type='node',
                        label=diagram[x].name,
                        degenerate=dim == diagram[x].dim)
        if dim >= 1:
            for x in diagram.shape[dim-1]:
                base.add_node(
                        x,
                        type='wire',
                        label=diagram[x].name,
                        degenerate=dim-1 == diagram[x].dim)

        vert = base.copy()
        if dim >= 1:
            for x in diagram.shape[dim]:
                for y in diagram.shape.faces(x, '-'):
                    vert.add_edge(y, x)
                for y in diagram.shape.faces(x, '+'):
                    vert.add_edge(x, y)
        horiz = base.copy()
        if dim >= 2:
            for x in horiz.nodes:
                for y in horiz.nodes:
                    if y != x:
                        out_x = GrSubset(
                            GrSet(x), diagram.shape,
                            wfcheck=False).closure().boundary_max(
                                    '+', dim-2)
                        in_y = GrSubset(
                            GrSet(y), diagram.shape,
                            wfcheck=False).closure().boundary_max(
                                    '-', dim-2)
                        if not out_x.isdisjoint(in_y):
                            horiz.add_edge(x, y)
        depth = base.copy()
        if dim >= 3:
            for x in diagram.shape[dim-1]:
                for y in diagram.shape[dim-1]:
                    if y != x:
                        out_x = GrSubset(
                            GrSet(x), diagram.shape,
                            wfcheck=False).closure().boundary_max(
                                    '+', dim-3)
                        in_y = GrSubset(
                            GrSet(y), diagram.shape,
                            wfcheck=False).closure().boundary_max(
                                    '-', dim-3)
                        if not out_x.isdisjoint(in_y):
                            depth.add_edge(x, y)

        def remove_cycles(graph):
            cycles = list(nx.simple_cycles(graph))
            to_delete = set()
            for cycle in cycles:
                for i in range(len(cycle) - 1):
                    to_delete.add((cycle[i], cycle[i+1]))
                to_delete.add((cycle[-1], cycle[0]))
            graph.remove_edges_from(to_delete)
        remove_cycles(horiz)
        remove_cycles(depth)

        self._vert = vert
        self._horiz = horiz
        self._depth = depth
