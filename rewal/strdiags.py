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
        if isinstance(diagram, Diagram):
            shape = diagram.shape
            generators = diagram.ambient.generators
            self.name = diagram.name

            def isdegenerate(x):
                return generators[diagram[x]]['shape'].dim != x.dim
        else:
            if isinstance(diagram, rewal.shapes.Shape):
                diagram = diagram.id()
            if isinstance(diagram, rewal.shapes.ShapeMap):
                shape = diagram.source
                generators = {x: {} for x in diagram.target}
                self.name = str(diagram)

                def isdegenerate(x):
                    return diagram[x].dim != x.dim
            else:
                raise TypeError(utils.type_err(
                    Diagram, diagram))

        dim = shape.dim

        self._nodes = {
                x: {
                    'label': diagram[x],
                    'color': generators[diagram[x]].get(
                        'color', None),
                    'stroke': generators[diagram[x]].get(
                        'stroke',
                        generators[diagram[x]].get(
                            'color', None)),
                    'draw_node': generators[diagram[x]].get(
                        'draw_node', True),
                    'draw_label': generators[diagram[x]].get(
                        'draw_label', True),
                    'isdegenerate': isdegenerate(x)
                    }
                for x in shape[dim]}
        self._wires = {
                x: {
                    'label': diagram[x],
                    'color': generators[diagram[x]].get(
                        'stroke',
                        generators[diagram[x]].get(
                            'color', None)),
                    'draw_label': generators[diagram[x]].get(
                        'draw_label', True),
                    'isdegenerate': isdegenerate(x)
                    }
                for x in shape[dim-1]}

        graph = nx.DiGraph()
        graph.add_nodes_from(self.nodes)
        graph.add_nodes_from(self.wires)
        if dim >= 1:
            out_1 = dict()
            in_1 = dict()
            for x in self.nodes:
                out_1[x] = shape.faces(x, '+')
                in_1[x] = shape.faces(x, '-')
                for y in in_1[x]:
                    graph.add_edge(y, x)
                for y in out_1[x]:
                    graph.add_edge(x, y)

        widthgraph = nx.DiGraph()
        widthgraph.add_nodes_from(self.nodes)
        widthgraph.add_nodes_from(self.wires)
        if dim >= 2:
            out_2 = dict()
            in_2 = dict()
            for x in self.wires:
                out_2[x] = shape.faces(x, '+')
                in_2[x] = shape.faces(x, '-')
            for x in self.nodes:
                out_2[x] = rewal.ogposets.GrSet(*[
                    z for w in out_1[x]
                    for z in shape.faces(w, '+')
                    if shape.cofaces(z, '-').isdisjoint(out_1[x])])
                in_2[x] = rewal.ogposets.GrSet(*[
                    z for w in in_1[x]
                    for z in shape.faces(w, '-')
                    if shape.cofaces(z, '+').isdisjoint(in_1[x])])
            for x in widthgraph:
                for y in widthgraph:
                    if not out_2[x].isdisjoint(in_2[y]):
                        widthgraph.add_edge(x, y)

        depthgraph = nx.DiGraph()
        depthgraph.add_nodes_from(self.wires)
        if dim >= 3:
            out_3 = dict()
            in_3 = dict()
            for x in depthgraph:
                out_3[x] = rewal.ogposets.GrSet(*[
                    z for w in out_2[x]
                    for z in shape.faces(w, '+')
                    if shape.cofaces(z, '-').isdisjoint(out_2[x])])
                in_3[x] = rewal.ogposets.GrSet(*[
                    z for w in in_2[x]
                    for z in shape.faces(w, '-')
                    if shape.cofaces(z, '+').isdisjoint(in_2[x])])
            for x in depthgraph:
                for y in depthgraph:
                    if not out_3[x].isdisjoint(in_3[y]):
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

    def __str__(self):
        return '{} with {} nodes and {} wires'.format(
                type(self).__name__, str(len(self.nodes)),
                str(len(self.wires)))

    def __eq__(self, other):
        return isinstance(other, StrDiag) and \
                self.graph == other.graph and \
                self.widthgraph == other.widthgraph and \
                self.depthgraph == other.depthgraph

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
        sources, sinks = [], []
        for x in self.graph:
            coordinates[x] = (
                    (longest_width[x][0] + 1) / (sum(longest_width[x]) + 2),
                    (longest_height[x][0] + 1) / (sum(longest_height[x]) + 2)
                    )
            if self.graph.in_degree(x) == 0:
                sources.append(x)
            if self.graph.out_degree(x) == 0:
                sinks.append(x)

        def solve_clashes(coord_dict):
            for coord in set(coord_dict.values()):  # Solve clashes
                keys = [x for x in coord_dict if coord_dict[x] == coord]
                if len(keys) > 1:
                    n = len(keys)
                    for k, x in enumerate(keys):
                        coordinates[x] = (
                            coordinates[x][0] + (xstep/3)*cos(.4 + (2*pi*k)/n),
                            coordinates[x][1] + (ystep/3)*sin(.4 + (2*pi*k)/n)
                            )
        solve_clashes(coordinates)  # xy clashes in initial placement

        sources_x = {x: coordinates[x][0] for x in sources}
        sinks_x = {x: coordinates[x][0] for x in sinks}
        solve_clashes(sources_x)  # x clashes in input boundary
        solve_clashes(sinks_x)  # x clashes in output boundary

        return coordinates

    def draw(self, **params):
        """
        Draws the string diagram with a backend.
        """
        # Parameters
        show = params.get('show', True)

        bgcolor = params.get('bgcolor', 'white')
        fgcolor = params.get('fgcolor', 'black')
        infocolor = params.get('infocolor', 'red')
        wirecolor = params.get('wirecolor', fgcolor)
        nodecolor = params.get('nodecolor', fgcolor)
        nodestroke = params.get('nodestroke', nodecolor)
        degenalpha = params.get('degenalpha', 0.1)

        labels = params.get('labels', True)
        wirelabels = params.get('wirelabels', labels)
        nodelabels = params.get('nodelabels', labels)
        labeloffset = params.get('labeloffset', (4, 4))

        positions = params.get('positions', False)
        wirepositions = params.get('wirepositions', positions)
        nodepositions = params.get('nodepositions', positions)
        positionoffset = params.get('positionoffset', (4, -16))

        orientation = params.get('orientation', 'bt')

        coord = self.place_vertices()
        backend = MatBackend(
                bgcolor=bgcolor,
                fgcolor=fgcolor,
                orientation=orientation,
                degenalpha=degenalpha,
                name=self.name)

        wiresort = list(nx.topological_sort(self.depthgraph))
        for wire in reversed(wiresort):
            color = wirecolor if self.wires[wire]['color'] is None \
                else self.wires[wire]['color']

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

            if wirepositions:
                backend.draw_label(
                        str(wire.pos),
                        coord[wire],
                        positionoffset,
                        color=infocolor)

            if wirelabels and self.wires[wire]['draw_label']:
                backend.draw_label(
                        self.wires[wire]['label'],
                        coord[wire],
                        labeloffset)

        def is_drawn(node):
            if self.nodes[node]['isdegenerate']:
                return False
            return self.nodes[node]['draw_node']

        for node in self.nodes:
            stroke = nodestroke if self.nodes[node]['stroke'] is None \
                else self.nodes[node]['stroke']
            color = nodecolor if self.nodes[node]['color'] is None \
                else self.nodes[node]['color']

            if is_drawn(node):
                backend.draw_node(
                    coord[node],
                    color,
                    stroke)
                if nodelabels and self.nodes[node]['draw_label']:
                    backend.draw_label(
                        self.nodes[node]['label'],
                        coord[node], labeloffset)

            if nodepositions:
                backend.draw_label(
                        str(node.pos),
                        coord[node],
                        positionoffset,
                        color=infocolor)

        if show:
            backend.show()


class DrawBackend(ABC):
    def __init__(self, **params):
        self.bgcolor = params.get('bgcolor')
        self.fgcolor = params.get('fgcolor')
        self.orientation = params.get('orientation')
        self.degenalpha = params.get('degenalpha')
        self.name = params.get('name')

    def rotate(self, xy):
        if self.orientation == 'tb':
            return (xy[0], 1-xy[1])
        if self.orientation == 'lr':
            return (xy[1], 1-xy[0])
        if self.orientation == 'rl':
            return (1-xy[1], 1-xy[0])
        return xy


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
        y_offset = 0
        # y_offset = .2*(wire_xy[1] - node_xy[1])
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
                alpha=self.degenalpha if isdegenerate else 1,
                lw=1)
        self.axes.add_patch(p_contour)
        self.axes.add_patch(p_wire)

    def draw_label(self, label, xy, offset, **params):
        color = params.get('color', self.fgcolor)
        size = params.get('size', 'medium')
        weight = params.get('weight', 'regular')

        xy = self.rotate(xy)
        xytext = (xy[0] + offset[0], xy[1] + offset[1])
        self.axes.annotate(
                label,
                xy,
                xytext=xytext,
                textcoords='offset pixels',
                color=color,
                fontsize=size,
                fontweight=weight)

    def draw_node(self, xy, color, stroke):
        xy = self.rotate(xy)
        self.axes.scatter(
                xy[0],
                xy[1],
                s=40,
                c=color,
                edgecolors=stroke)

    def show(self):
        self.fig.subplots_adjust(
                top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
        self.fig.canvas.manager.set_window_title(self.name)
        self.fig.show()


def draw(*diagrams, **params):
    for diagram in diagrams:
        StrDiag(diagram).draw(**params)


def draw_boundaries(diagram, dim=None, **params):
    """
    Draws the input and the output boundaries of a diagram.
    """
    StrDiag(diagram.boundary('-', dim)).draw(**params)
    StrDiag(diagram.boundary('+', dim)).draw(**params)
