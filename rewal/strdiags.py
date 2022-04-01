"""
Implements string diagram visualisations.
"""

import networkx as nx
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch

from rewal import utils
from rewal.ogposets import (GrSet, GrSubset)
from rewal.diagrams import Diagram

matplotlib.use('TkCairo')


class StrDiag:
    """
    Class for string diagrams.
    """
    def __init__(self, diagram, **kwargs):
        utils.typecheck(diagram, {'type': Diagram})
        dim = diagram.dim

        self._nodes = {
                x: {
                    'label': diagram[x].name,
                    'color': diagram.ambient.generators[diagram[x].name].get(
                        'color', 'black'),
                    'isdegenerate': dim != diagram[x].dim
                    }
                for x in diagram.shape[dim]}
        self._wires = {
                x: {
                    'label': diagram[x].name,
                    'color': diagram.ambient.generators[diagram[x].name].get(
                        'color', 'black'),
                    'isdegenerate': dim-1 != diagram[x].dim
                    }
                for x in diagram.shape[dim-1]}

        graph = nx.DiGraph()
        graph.add_nodes_from(self.nodes)
        graph.add_nodes_from(self.wires)
        if dim >= 1:
            for x in self.nodes:
                for y in diagram.shape.faces(x, '-'):
                    graph.add_edge(y, x)
                for y in diagram.shape.faces(x, '+'):
                    graph.add_edge(x, y)

        width = nx.DiGraph()
        width.add_nodes_from(self.nodes)
        width.add_nodes_from(self.wires)
        if dim >= 2:
            for x in width:
                for y in width:
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
                            width.add_edge(x, y)

        depth = nx.DiGraph()
        depth.add_nodes_from(self.wires)
        if dim >= 3:
            for x in depth:
                for y in depth:
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

        remove_cycles(width)
        remove_cycles(depth)

        self._graph = graph
        self._width = width
        self._depth = depth

    @property
    def graph(self):
        return self._graph

    @property
    def width(self):
        return self._width

    @property
    def depth(self):
        return self._depth

    @property
    def nodes(self):
        return self._nodes

    @property
    def wires(self):
        return self._wires

    def place_vertices(self):
        """ Places vertices on unit square. """
        def longest_paths(graph):
            tsort = list(nx.topological_sort(graph))
            longest_paths = dict()
            for x in tsort:
                longest_fw = {y: -1 for y in graph}
                longest_fw[x] = 0
                for y in (y for y in tsort if longest_fw[y] >= 0):
                    for z in graph.successors(y):
                        if longest_fw[z] < longest_fw[y] + 1:
                            longest_fw[z] = longest_fw[y] + 1
                longest_bw = {y: -1 for y in graph}
                longest_bw[x] = 0
                for y in (y for y in reversed(tsort) if longest_bw[y] >= 0):
                    for z in graph.predecessors(y):
                        if longest_bw[z] < longest_bw[y] + 1:
                            longest_bw[z] = longest_bw[y] + 1
                longest_paths[x] = (
                        max(longest_bw.values()),
                        max(longest_fw.values()))
            return longest_paths
        longest_height = longest_paths(self.graph)
        longest_width = longest_paths(self.width)

        coordinates = dict()
        for x in self.graph:
            coordinates[x] = (
                    (longest_width[x][0] + 0.5) / (sum(longest_width[x]) + 1),
                    (longest_height[x][0] + 0.5) / (sum(longest_height[x]) + 1)
                    )
        return coordinates

    def draw(self):  # Just a stub to see if all works.
        ax = plt.subplots()[1]
        coord = self.place_vertices()
        for x in self.wires:
            for y in [
                    *self.graph.predecessors(x),
                    *self.graph.successors(x)
                    ]:
                path = Path(
                        [
                            coord[x],
                            (coord[x][0], coord[y][1]),
                            coord[y]
                        ], [
                            Path.MOVETO,
                            Path.CURVE3,
                            Path.CURVE3
                        ])
                patch = PathPatch(
                        path,
                        facecolor='none',
                        edgecolor=self.wires[x]['color'],
                        alpha=0.2 if self.wires[x]['isdegenerate'] else 1,
                        lw=1)
                ax.add_patch(patch)
            if self.graph.in_degree(x) == 0:
                path = Path(
                        [coord[x], (coord[x][0], 0)],
                        [Path.MOVETO, Path.LINETO])
                patch = PathPatch(
                        path,
                        facecolor='none',
                        edgecolor=self.wires[x]['color'],
                        alpha=0.2 if self.wires[x]['isdegenerate'] else 1,
                        lw=1)
                ax.add_patch(patch)
            if self.graph.out_degree(x) == 0:
                path = Path(
                        [coord[x], (coord[x][0], 1)],
                        [Path.MOVETO, Path.LINETO])
                patch = PathPatch(
                        path,
                        facecolor='none',
                        edgecolor=self.wires[x]['color'],
                        alpha=0.2 if self.wires[x]['isdegenerate'] else 1,
                        lw=1)
                ax.add_patch(patch)

        for x in (
                x for x in self.nodes
                if not self.nodes[x]['isdegenerate']):
            ax.scatter(
                    coord[x][0], coord[x][1],
                    c=self.nodes[x]['color'])

        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        plt.axis('off')
        plt.subplots_adjust(
                top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
        plt.show()
