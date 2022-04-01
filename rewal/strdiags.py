"""
Implements string diagram visualisations.
"""

import networkx as nx

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patheffects as patheffects
from matplotlib.path import Path
from matplotlib.patches import PathPatch

from math import (sin, cos, pi)

import rewal
from rewal import utils
from rewal.diagrams import Diagram

matplotlib.use('TkCairo')


class StrDiag:
    """
    Class for string diagrams.
    """
    def __init__(self, diagram, **kwargs):
        if isinstance(diagram, rewal.shapes.ShapeMap):
            diagram = Diagram.yoneda(diagram)
        else:
            if isinstance(diagram, rewal.shapes.Shape):
                diagram = Diagram.yoneda(diagram.id())
            else:
                utils.typecheck(diagram, {'type': Diagram})

        dim = diagram.dim
        generators = diagram.ambient.generators

        self._nodes = {
                x: {
                    'label': diagram[x].name,
                    'color': generators[diagram[x].name].get(
                        'color', 'black'),
                    'draw_node': generators[diagram[x].name].get(
                        'draw_node', True),
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

        widthgraph = nx.DiGraph()
        widthgraph.add_nodes_from(self.nodes)
        widthgraph.add_nodes_from(self.wires)
        if dim >= 2:
            for x in widthgraph:
                for y in widthgraph:
                    if y != x:
                        out_x = rewal.ogposets.GrSubset(
                            rewal.ogposets.GrSet(x), diagram.shape,
                            wfcheck=False).closure().boundary_max(
                                    '+', dim-2)
                        in_y = rewal.ogposets.GrSubset(
                            rewal.ogposets.GrSet(y), diagram.shape,
                            wfcheck=False).closure().boundary_max(
                                    '-', dim-2)
                        if not out_x.isdisjoint(in_y):
                            widthgraph.add_edge(x, y)

        depthgraph = nx.DiGraph()
        depthgraph.add_nodes_from(self.wires)
        if dim >= 3:
            for x in depthgraph:
                for y in depthgraph:
                    if y != x:
                        out_x = rewal.ogposets.GrSubset(
                            rewal.ogposets.GrSet(x), diagram.shape,
                            wfcheck=False).closure().boundary_max(
                                    '+', dim-3)
                        in_y = rewal.ogposets.GrSubset(
                            rewal.ogposets.GrSet(y), diagram.shape,
                            wfcheck=False).closure().boundary_max(
                                    '-', dim-3)
                        if not out_x.isdisjoint(in_y):
                            depthgraph.add_edge(x, y)

        def remove_cycles(graph):
            cycles = list(nx.simple_cycles(graph))
            to_delete = set()
            for cycle in cycles:
                for i in range(len(cycle) - 1):
                    to_delete.add((cycle[i], cycle[i+1]))
                to_delete.add((cycle[-1], cycle[0]))
            graph.remove_edges_from(to_delete)

        remove_cycles(widthgraph)
        remove_cycles(depthgraph)

        self._graph = graph
        self._widthgraph = widthgraph
        self._depthgraph = depthgraph

    @property
    def graph(self):
        return self._graph

    @property
    def widthgraph(self):
        return self._widthgraph

    @property
    def depthgraph(self):
        return self._depthgraph

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
        longest_width = longest_paths(self.widthgraph)
        longest_height = longest_paths(self.graph)

        xstep = 1 / (max(
            [longest_width[x][0] for x in longest_width]) + 2)
        ystep = 1 / (max(
            [longest_height[x][0] for x in longest_height]) + 2)

        coordinates = dict()
        for x in self.graph:
            coordinates[x] = (
                    (longest_width[x][0] + 1) / (sum(longest_width[x]) + 2),
                    (longest_height[x][0] + 1) / (sum(longest_height[x]) + 2)
                    )

        for coord in set(coordinates.values()):  # Solve clashes
            keys = [x for x in coordinates if coordinates[x] == coord]
            if len(keys) > 1:
                n = len(keys)
                for k, x in enumerate(keys):
                    coordinates[x] = (
                        coordinates[x][0] + (xstep/4)*cos(.4 + (2*pi*k)/n),
                        coordinates[x][1] + (ystep/4)*sin(.4 + (2*pi*k)/n)
                        )
        return coordinates

    def draw(self):  # Just a stub to see if all works.
        ax = plt.subplots()[1]
        coord = self.place_vertices()

        wiresort = list(nx.topological_sort(self.depthgraph))

        contour = [
                patheffects.Stroke(linewidth=3, alpha=1, foreground='white'),
                patheffects.Normal()
                ]
        for x in reversed(wiresort):
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
                        alpha=0.1 if self.wires[x]['isdegenerate'] else 1,
                        path_effects=contour,
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
                        alpha=0.1 if self.wires[x]['isdegenerate'] else 1,
                        path_effects=contour,
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
                        alpha=0.1 if self.wires[x]['isdegenerate'] else 1,
                        path_effects=contour,
                        lw=1)
                ax.add_patch(patch)
            ax.annotate(
                    self.wires[x]['label'],
                    (coord[x][0]+.01, coord[x][1])
                    )

        def is_drawn(x):
            if self.nodes[x]['isdegenerate']:
                return False
            return self.nodes[x]['draw_node']

        for x in (
                x for x in self.nodes
                if is_drawn(x)):
            ax.scatter(
                    coord[x][0], coord[x][1],
                    s=40,
                    c=self.nodes[x]['color'])
            ax.annotate(
                    self.nodes[x]['label'],
                    (coord[x][0]+.01, coord[x][1]+.01)
                    )

        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        plt.axis('off')
        plt.subplots_adjust(
                top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
        plt.show()
