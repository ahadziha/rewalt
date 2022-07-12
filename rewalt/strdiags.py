"""
Implements string diagram visualisations.
"""

import networkx as nx

from rewalt import (utils, ogposets, shapes, diagrams, drawing)


DEFAULT = {
        'tikz': False,
        'scale': 3,
        'show': True,
        'depth': True,
        'bgcolor': 'white',
        'fgcolor': 'black',
        'infocolor': 'magenta',
        'degenalpha': 0.1,
        'labels': True,
        'labeloffset': (4, 4),
        'positions': False,
        'positionoffset': (4, -16),
        'positionoffsettikz': (4, -6),
        'orientation': 'bt'}


class StrDiag:
    """
    Class for string diagram visualisations of diagrams and shapes.

    A string diagram depicts a top-dimensional "slice" of a diagram.
    The top-dimensional cells are represented as *nodes*, and the
    codimension-1 cells are represented as *wires*. The inputs of a
    top-dimensional cell are incoming wires of the associated node,
    and the outputs are outgoing wires.

    The input->node->output order determines an acyclic flow
    between nodes and wires, which is represented in a string diagram
    by placing them at different "heights".

    There are two other "flows" that we take into account:

    - from codimension-2 inputs, to top-dimensional or codimension-1
      cell, to codimension-2 outputs (only in dimension > 1);
    - from codimension-3 inputs, to codimension-1 cells, to
      codimension-3 outputs (only in dimension > 2).

    These are not in general acyclic; however, we obtain an acyclic
    flow by removing all directed loops. If there is a flow of the first
    kind between nodes and wires, we place them at different "widths".

    If there is a flow of the second kind between wires, we place them
    at different "depths"; this is only seen when wires cross each other,
    in which case the one of lower depth is depicted as passing over
    the one of higher depth.

    Internally, these data are encoded as a triple of NetworkX directed
    graphs, sharing the same vertices, partitioned into "node vertices"
    and "wire vertices". These graphs encode the "main (height) flow", the
    "width flow" and the "depth flow" between nodes and wires.

    The class then contains a method :meth:`place_vertices` that places
    the vertices on a [0, 1]x[0, 1] canvas, taking into account the
    height and width relations and resolving clashes.

    Finally, it contains a method :meth:`draw` that outputs a
    visualisation of the string diagram. The visualisation has
    customisable colours, orientation, and labels, and works with any
    :class:`drawing.DrawBackend`; currently available are

    - a Matplotlib backend, and
    - a TikZ backend.

    Arguments
    ---------
    diagram : :class:`diagrams.Diagram | shapes.Shape | shapes.ShapeMap`
        A diagram or a shape or a shape map.

    Notes
    -----
    The "main flow" graph is essentially the *open graph* encoding of
    the string diagram in the sense of Dixon & Kissinger.
    """
    def __init__(self, diagram):
        if isinstance(diagram, diagrams.Diagram):
            shape = diagram.shape
            generators = diagram.ambient.generators
            self.name = diagram.name

            def isdegenerate(x):
                return generators[diagram[x]]['shape'].dim != x.dim
        else:
            if isinstance(diagram, shapes.Shape):
                diagram = diagram.id()
            if isinstance(diagram, shapes.ShapeMap):
                shape = diagram.source
                generators = {x: {} for x in diagram.target}
                self.name = str(diagram)

                def isdegenerate(x):
                    return diagram[x].dim != x.dim
            else:
                raise TypeError(utils.type_err(
                    diagrams.Diagram, diagram))

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
                out_2[x] = ogposets.GrSet(*[
                    z for w in out_1[x]
                    for z in shape.faces(w, '+')
                    if shape.cofaces(z, '-').isdisjoint(out_1[x])])
                in_2[x] = ogposets.GrSet(*[
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
                out_3[x] = ogposets.GrSet(*[
                    z for w in out_2[x]
                    for z in shape.faces(w, '+')
                    if shape.cofaces(z, '-').isdisjoint(out_2[x])])
                in_3[x] = ogposets.GrSet(*[
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
        """
        Returns the main flow graph between node and wire vertices.

        Returns
        -------
        graph : :class:`networkx.DiGraph`
            The main flow graph.
        """
        return self._graph

    @property
    def widthgraph(self):
        """
        Returns the "width" flow graph between node and wire vertices.

        Returns
        -------
        widthgraph : :class:`networkx.DiGraph`
            The width flow graph.
        """
        return self._widthgraph

    @property
    def depthgraph(self):
        """
        Returns the "depth" flow graph between wire vertices.

        Returns
        -------
        depthgraph : :class:`networkx.DiGraph`
            The depth flow graph.
        """
        return self._depthgraph

    @property
    def nodes(self):
        """
        Returns the nodes of the string diagram, together with all
        the stored associated information.

        This is a dictionary whose keys are the elements
        of the diagram's shape corresponding to nodes. For each node, the
        object stores another dictionary, which contains

        - the node's label (:code:`label`),
        - the node's fill colour (:code:`color`) and stroke colour
          (:code:`stroke`),
        - booleans specifying whether to draw the node and/or its label
          (:code:`draw_node`, :code:`draw_label`), and
        - a boolean specifying whether the node represents a degenerate
          cell (:code:`isdegenerate`).

        Returns
        -------
        nodes : :class:`dict[dict]`
            The nodes of the string diagram.
        """
        return self._nodes

    @property
    def wires(self):
        """
        Returns the wires of the string diagram, together with all
        the stored associated information.

        This is a dictionary whose keys are the elements
        of the diagram's shape corresponding to wires. For each node, the
        object stores another dictionary, which contains

        - the wire's label (:code:`label`),
        - the wire's colour (:code:`color`),
        - a boolean specifying whether to draw the wire's label
          (:code:`draw_label`), and
        - a boolean specifying whether the wire represents a degenerate
          cell (:code:`isdegenerate`).

        Returns
        -------
        wires : :class:`dict[dict]`
            The nodes of the string diagram.
        """
        return self._wires

    def place_vertices(self):
        """
        Places node and wire vertices on the unit square canvas, and
        returns their coordinates.

        The node and wire vertices are first placed on different heights
        and widths, proportional to the ratio between the longest path
        to the vertex and the longest path from the vertex in the main
        flow graph and the width flow graph.

        In dimension > 2, this may result in clashes, where some vertices
        are given the same coordinates. In this case, these are
        resolved by "splitting" the clashing vertices, placing them
        at equally spaced angles of a circle centred on the clash
        coordinates, with an appropriately small radius that does not
        result in further clashes.

        The coordinates are returned as a dictionary whose keys are
        the elements corresponding to nodes and wires.

        Returns
        -------
        coordinates : :class:`dict[tuple[float]]`
            The coordinates assigned to wire and node vertices.
        """
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
            [longest_width[x][0] for x in longest_width], default=0) + 2)
        ystep = 1 / (max(
            [longest_height[x][0] for x in longest_height], default=0) + 2)

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
            from math import (sin, cos, pi)

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
        Outputs a visualisation of the string diagram, using a backend.

        Currently supported are a Matplotlib backend and a TikZ backend;
        in both cases it is possible to show the output (as a pop-up
        window for Matplotlib, or as code for TikZ) or save to file.

        Various customisation options are available, including different
        orientations and colours.

        Keyword arguments
        -----------------
        tikz : :class:`bool`
            Whether to output TikZ code (default is :code:`False`).
        show : :class:`bool`
            Whether to show the output (default is :code:`True`).
        path : :class:`str`
            Path where to save the output (default is :code:`None`).
        orientation : :class:`str`
            Orientation of the string diagram: one of :code:`'bt'`
            (bottom-to-top), :code:`'lr'` (left-to-right),
            :code:`'tb'` (top-to-bottom), :code:`'rl'` (right-to-left)
            (default is :code:`'bt'`).
        depth : :class:`bool`
            Whether to take into account the depth flow graph when
            drawing wires (default is :code:`True`).
        bgcolor : multiple types
            The background colour (default is :code:`'white'`).
        fgcolor : multiple types
            The foreground colour, given by default to nodes, wires,
            and labels (default is :code:`'black'`).
        infocolor : multiple types
            The colour of additional information displayed in
            the diagram, such as positions (default is :code:`'magenta'`).
        wirecolor : multiple types
            The default wire colour (default is same as `fgcolor`).
        nodecolor : multiple types
            The default node fill colour (default is same as `fgcolor`).
        nodestroke : multiple types
            The default node stroke colour (default is same as `nodecolor`).
        degenalpha : :class:`float`
            The alpha factor of wires corresponding to degenerate cells
            (default is :code:`0.1`).
        labels : :class:`bool`
            Whether to display node and wire labels (default is
            :code:`True`).
        nodelabels : :class:`bool`
            Whether to display node labels (default is same as `labels`).
        wirelabels : :class:`bool`
            Whether to display wire labels (default is same as `labels`).
        labeloffset : :class:`tuple[float]`
            Point offset of labels relative to vertices (default is
            :code:`(4, 4)`).
        positions : :class:`bool`
            Whether to display node and wire positions (default is
            :code:`False`).
        nodepositions : :class:`bool`
            Whether to display node positions (default is same as
            `positions`).
        wirepositions : :class:`bool`
            Whether to display wire positions (default is same as
            `positions`).
        positionoffset : :class:`tuple[float]`
            Point offset of positions relative to vertices (default is
            :code:`(4, -16)` for Matplotlib, :code:`(4, -6)` for TikZ).
        scale : :class:`float`
            (TikZ only) Scale factor to apply to output (default is
            :code:`3`).
        xscale : :class:`float`
            (TikZ only) Scale factor to apply to x axis in output
            (default is same as `scale`)
        yscale : :class:`float`
            (TikZ only) Scale factor to apply to y axis in output
            (default is same as `scale`)
        """
        # Parameters
        tikz = params.get('tikz', DEFAULT['tikz'])
        scale = params.get('scale', DEFAULT['scale'])
        xscale = params.get('xscale', scale)
        yscale = params.get('yscale', scale)

        show = params.get('show', DEFAULT['show'])
        path = params.get('path', None)

        depth = params.get('depth', DEFAULT['depth'])
        bgcolor = params.get(
                'bgcolor', DEFAULT['bgcolor'])
        fgcolor = params.get(
                'fgcolor', DEFAULT['fgcolor'])
        infocolor = params.get(
                'infocolor', DEFAULT['infocolor'])
        wirecolor = params.get('wirecolor', fgcolor)
        nodecolor = params.get('nodecolor', fgcolor)
        nodestroke = params.get('nodestroke', nodecolor)
        degenalpha = params.get(
                'degenalpha', DEFAULT['degenalpha'])

        labels = params.get(
                'labels', DEFAULT['labels'])
        wirelabels = params.get('wirelabels', labels)
        nodelabels = params.get('nodelabels', labels)
        labeloffset = params.get(
                'labeloffset', DEFAULT['labeloffset'])

        positions = params.get(
                'positions', DEFAULT['positions'])
        wirepositions = params.get('wirepositions', positions)
        nodepositions = params.get('nodepositions', positions)

        podefault = DEFAULT['positionoffsettikz'] \
            if tikz else DEFAULT['positionoffset']
        positionoffset = params.get(
                'positionoffset', podefault)

        orientation = params.get(
                'orientation', DEFAULT['orientation'])

        coord = self.place_vertices()

        backendclass = drawing.TikZBackend if tikz else drawing.MatBackend
        backend = backendclass(
                bgcolor=bgcolor,
                fgcolor=fgcolor,
                orientation=orientation,
                name=self.name)

        wiresort = list(nx.topological_sort(self.depthgraph))
        for wire in reversed(wiresort):
            color = wirecolor if self.wires[wire]['color'] is None \
                else self.wires[wire]['color']
            alpha = degenalpha if self.wires[wire]['isdegenerate'] else 1

            for node in [
                    *self.graph.predecessors(wire),
                    *self.graph.successors(wire)
                    ]:
                backend.draw_wire(
                        coord[wire], coord[node],
                        color=color,
                        alpha=alpha,
                        depth=depth)

            if self.graph.in_degree(wire) == 0:
                backend.draw_wire(
                        coord[wire], (coord[wire][0], 0),
                        color=color,
                        alpha=alpha,
                        depth=depth)

            if self.graph.out_degree(wire) == 0:
                backend.draw_wire(
                        coord[wire], (coord[wire][0], 1),
                        color=color,
                        alpha=alpha,
                        depth=depth)

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
                    color=color,
                    stroke=stroke)
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

        backend.output(
                path=path,
                show=show,
                xscale=xscale,
                yscale=yscale)


def draw(*diagrams, **params):
    """
    Given any number of diagrams, generates their string
    diagrams and draws them.

    This is the same as generating the string diagram for each
    diagram, and calling :meth:`StrDiag.draw` with the given
    parameters on each one of them.

    Arguments
    ---------
    *diagrams : :class:`diagrams.Diagram | shapes.Shape | shapes.ShapeMap`
        Any number of diagrams or shapes or shape maps.

    Keyword arguments
    -----------------
    **params
        Passed to :meth:`StrDiag.draw`.
    """
    for diagram in diagrams:
        StrDiag(diagram).draw(**params)


def draw_boundaries(diagram, dim=None, **params):
    """
    Given a diagram, generates the string diagram of its input and
    output boundaries of a given dimension, and draws them.

    Arguments
    ---------
    diagram : :class:`diagrams.Diagram | shapes.Shape | shapes.ShapeMap`
        A diagram or a shape or a shape map.
    dim : :class:`int`, optional
        Dimension of the boundary (default is :code:`diagram.dim - 1`).

    Keyword arguments
    -----------------
    *params
        Passed to :meth:`StrDiag.draw`.
    """
    StrDiag(diagram.boundary('-', dim)).draw(**params)
    StrDiag(diagram.boundary('+', dim)).draw(**params)


def to_gif(diagram, *diagrams, **params):
    """
    Given a non-zero number of diagrams, generates their string
    diagrams and outputs a GIF animation of the sequence of their
    visualisations.

    Arguments
    ---------
    diagram : :class:`diagrams.Diagram | shapes.Shape | shapes.ShapeMap`
        A diagram or a shape or a shape map.
    *diagrams : :class:`diagrams.Diagram | shapes.Shape | shapes.ShapeMap`
        Any number of diagrams or shapes or shape maps.

    Keyword arguments
    -----------------
    timestep : :class:`int`
        The time step for the animation (default is :code:`1000`).
    loop : :class:`bool`
        Whether to loop around the animation (default is :code:`True`).
    **params
        Passed to :meth:`StrDiag.draw`.
    """
    import os
    from tempfile import NamedTemporaryFile, TemporaryDirectory
    from PIL import Image

    path = params.pop('path', None)
    params.pop('show', False)
    timestep = params.get('timestep', 1000)
    loop = params.get('loop', False)
    frames = []

    path = path or os.path.basename(NamedTemporaryFile(
        suffix='.gif', prefix='tmp_', dir='.').name)
    with TemporaryDirectory() as directory:
        for k, step in enumerate((diagram, *diagrams)):
            tmp_path = os.path.join(directory, '{}.png'.format(k))
            draw(step, path=tmp_path, show=False, **params)
            frames.append(Image.open(tmp_path))
        if loop:
            frames = frames + frames[::-1]
        frames[0].save(
                path,
                format='GIF',
                append_images=frames[1:],
                save_all=True,
                duration=timestep,
                **{'loop': 0} if loop else {})
