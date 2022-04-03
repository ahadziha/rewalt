"""
Implements string diagram visualisations.
"""

from abc import ABC

import networkx as nx

import matplotlib
import matplotlib.pyplot as plt
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
    def __init__(self, diagram, **params):
        if isinstance(diagram, rewal.shapes.ShapeMap):
            diagram = Diagram.yoneda(diagram)
        else:
            if isinstance(diagram, rewal.shapes.Shape):
                diagram = Diagram.yoneda(diagram.id())
            else:
                utils.typecheck(diagram, {'type': Diagram})

        dim = diagram.dim
        generators = diagram.ambient.generators

        self.name = diagram.name

        self._nodes = {
                x: {
                    'label': diagram[x].name,
                    'color': generators[diagram[x].name].get(
                        'color', None),
                    'draw_node': generators[diagram[x].name].get(
                        'draw_node', True),
                    'isdegenerate': dim != diagram[x].dim
                    }
                for x in diagram.shape[dim]}
        self._wires = {
                x: {
                    'label': diagram[x].name,
                    'color': generators[diagram[x].name].get(
                        'color', None),
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
                        coordinates[x][0] + (xstep/3)*cos(.4 + (2*pi*k)/n),
                        coordinates[x][1] + (ystep/3)*sin(.4 + (2*pi*k)/n)
                        )
        return coordinates

    def draw(self, **params):
        """
        Draws the string diagram with a backend.
        """
        # Parameters
        show = params.get('show', True)
        fgcolor = params.get('fgcolor', 'black')
        bgcolor = params.get('bgcolor', 'white')
        labels = params.get('labels', True)
        wire_labels = params.get('wire_labels', labels)
        node_labels = params.get('node_labels', labels)
        orientation = params.get('orientation', 'bt')

        coord = self.place_vertices()
        backend = MatBackend(
                bgcolor=bgcolor,
                fgcolor=fgcolor,
                orientation=orientation,
                name=self.name)

        wiresort = list(nx.topological_sort(self.depthgraph))
        for wire in reversed(wiresort):
            if self.wires[wire]['color'] is None:
                color = fgcolor
            else:
                color = self.wires[wire]['color']
            for node in [
                    *self.graph.predecessors(wire),
                    *self.graph.successors(wire)
                    ]:
                backend.draw_wire(
                        coord[wire], coord[node],
                        color,
                        self.wires[wire]['isdegenerate'])

            if self.graph.in_degree(wire) == 0:
                backend.draw_wire(
                        coord[wire], (coord[wire][0], 0),
                        color,
                        self.wires[wire]['isdegenerate'])

            if self.graph.out_degree(wire) == 0:
                backend.draw_wire(
                        coord[wire], (coord[wire][0], 1),
                        color,
                        self.wires[wire]['isdegenerate'])

            if wire_labels:
                backend.draw_label(
                        self.wires[wire]['label'],
                        coord[wire], (.01, .01))

        def is_drawn(node):
            if self.nodes[node]['isdegenerate']:
                return False
            return self.nodes[node]['draw_node']

        for node in (
                node for node in self.nodes
                if is_drawn(node)):
            if self.nodes[node]['color'] is None:
                color = fgcolor
            else:
                color = self.nodes[node]['color']
            backend.draw_node(
                    coord[node],
                    color)

            if node_labels:
                backend.draw_label(
                        self.nodes[node]['label'],
                        coord[node], (.01, .01))
        if show:
            backend.show()


class DrawBackend(ABC):
    def __init__(self, **params):
        self.bgcolor = params.get('bgcolor')
        self.fgcolor = params.get('fgcolor')
        self.orientation = params.get('orientation')
        self.name = params.get('name')


class MatBackend(DrawBackend):
    """
    Matplotlib drawing backend.
    """
    def __init__(self, **params):
        super().__init__(**params)

        self.fig, self.axes = plt.subplots()
        self.axes.set_facecolor(self.bgcolor)
        self.axes.set_xlim(0, 1)
        self.axes.set_ylim(0, 1)
        for side in ('top', 'right', 'bottom', 'left'):
            self.axes.spines[side].set_visible(False)

    def draw_wire(self, wire_xy, node_xy,
                  color, isdegenerate):
        """
        Draws a wire from a wire vertex to a node vertex.
        """
        y_offset = .2*(wire_xy[1] - node_xy[1])
        width = .02

        contour = Path(
                [
                    self.rotate(node_xy),
                    self.rotate(
                        (wire_xy[0] - (width/2), node_xy[1] + y_offset)),
                    self.rotate(
                        (wire_xy[0] - (width/2), wire_xy[1])),
                    self.rotate(
                        (wire_xy[0] + (width/2), wire_xy[1])),
                    self.rotate(
                        (wire_xy[0] + (width/2), node_xy[1] + y_offset)),
                    self.rotate(node_xy)
                ], [
                    Path.MOVETO,
                    Path.CURVE3,
                    Path.CURVE3,
                    Path.LINETO,
                    Path.CURVE3,
                    Path.CURVE3
                ])
        wire = Path(
                [
                    self.rotate(wire_xy),
                    self.rotate(
                        (wire_xy[0], node_xy[1] + y_offset)),
                    self.rotate(node_xy)
                ], [
                    Path.MOVETO,
                    Path.CURVE3,
                    Path.CURVE3
                ])
        p_contour = PathPatch(
                contour,
                facecolor=self.bgcolor,
                edgecolor='none')
        p_wire = PathPatch(
                wire,
                facecolor='none',
                edgecolor=color,
                alpha=0.1 if isdegenerate else 1,
                lw=1)
        self.axes.add_patch(p_contour)
        self.axes.add_patch(p_wire)

    def draw_label(self, label, xy, offset):
        xy = self.rotate(xy)
        xy = (xy[0] + offset[0], xy[1] + offset[1])
        self.axes.annotate(
                label,
                xy,
                color=self.fgcolor)

    def draw_node(self, xy, color):
        xy = self.rotate(xy)
        self.axes.scatter(
                xy[0],
                xy[1],
                s=40,
                c=color)

    def show(self):
        self.fig.subplots_adjust(
                top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
        self.fig.canvas.manager.set_window_title(self.name)
        self.fig.show()

    def rotate(self, xy):
        if self.orientation == 'tb':
            return (xy[0], 1-xy[1])
        if self.orientation == 'lr':
            return (xy[1], 1-xy[0])
        if self.orientation == 'rl':
            return (1-xy[1], 1-xy[0])
        return xy
