"""
Implements oriented Hasse diagram visualisation.
"""

import networkx as nx

from rewalt import utils, ogposets, diagrams, drawing


DEFAULT = {
        'tikz': False,
        'show': True,
        'bgcolor': 'white',
        'fgcolor': 'black',
        'labels': True,
        'inputcolor': 'magenta',
        'outputcolor': 'blue',
        'orientation': 'bt'}


class Hasse:
    """
    Class for "oriented Hasse diagrams" of oriented graded posets.

    The oriented Hasse diagram is stored as a NetworkX directed graph
    whose nodes are the elements of the oriented graded poset.

    The orientation information is encoded by having edges corresponding
    to input faces point *from* the face, and edges corresponding to
    output faces point *towards* the face. To recover the underlying
    poset's Hasse diagram, it suffices to reverse the edges that point
    from an element of higher dimension.

    Objects of the class can also store labels for nodes of the Hasse
    diagram, for example the images of the corresponding elements
    through a map or a diagram.

    The class also has a method :meth:`draw` that outputs a visualisation
    of the Hasse diagram. This works with any :class:`drawing.DrawBackend`;
    currently available are

    - a Matplotlib backend, and
    - a TikZ backend.

    Arguments
    ---------
    ogp : :class:`ogposets.OgPoset | ogposets.OgMap | diagrams.Diagram`
        The oriented graded poset, or a map of oriented graded posets,
        or a diagram.

    Notes
    -----
    If given a map of oriented graded posets (or shapes), produces the
    Hasse diagram of its source, with nodes labelled with the images
    of elements through the map.

    If given a diagram, produces the Hasse diagram of its shape, with
    nodes labelled with the images of elements through the diagram.
    """
    def __init__(self, ogp):
        if isinstance(ogp, ogposets.OgPoset):
            self._labels = ogp.id().mapping
        elif isinstance(ogp, diagrams.Diagram):
            self._labels = ogp.mapping
            ogp = ogp.shape
        elif isinstance(ogp, ogposets.OgMap):
            self._labels = ogp.mapping
            ogp = ogp.source
        else:
            raise TypeError(utils.type_err(
                ogposets.OgPoset, ogp))

        self._nodes = ogp.all().support

        diagram = nx.DiGraph()
        diagram.add_nodes_from(self.nodes)
        for x in self.nodes[1:]:
            for y in ogp.faces(x, '-'):
                diagram.add_edge(y, x, sign='-')
            for y in ogp.faces(x, '+'):
                diagram.add_edge(x, y, sign='+')

        self._diagram = diagram

    @property
    def nodes(self):
        """
        Returns the set of nodes of the Hasse diagram, that is, the
        graded set of elements of the oriented graded poset it encodes.

        Returns
        -------
        nodes : :class:`ogposets.GrSet`
            The set of nodes of the Hasse diagram.
        """
        return self._nodes

    @property
    def diagram(self):
        """
        Returns the oriented Hasse diagram as a NetworkX graph.

        Returns
        -------
        diagram : :class:`networkx.DiGraph`
            The oriented Hasse diagram.
        """
        return self._diagram

    @property
    def labels(self):
        """
        Returns the labels of nodes of the Hasse diagram, in the
        same format as :meth:`ogposets.OgMap.mapping`.

        Returns
        -------
        labels : :class:`list[list]`
            The labels of the Hasse diagram.
        """
        return self._labels

    def place_nodes(self):
        """
        Places the nodes of the Hasse diagram on the unit square
        canvas, and returns their coordinates.

        The nodes are placed on different heights according to the
        dimension of the element their correspond to.
        Elements of the same dimension are then placed at different
        widths in order of position.

        The coordinates are returned as a dictionary whose keys are
        the elements corresponding to nodes of the diagram.

        Returns
        -------
        coordinates : :class:`dict[tuple[float]]`
            The coordinates assigned to nodes.
        """
        dim = self.nodes.dim
        if dim < 0:
            return dict()

        ystep = 1/(dim+1)
        xstep = [1/(len(self.nodes[n])) for n in range(dim+1)]

        coordinates = dict()
        for x in self.nodes:
            coordinates[x] = (
                    (x.pos + 0.5)*xstep[x.dim],
                    (x.dim + 0.5)*ystep
                    )
        return coordinates

    def draw(self, **params):
        """
        Outputs a visualisation of the Hasse diagram, using a backend.

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
            Orientation of the Hasse diagram: one of :code:`'bt'`
            (bottom-to-top), :code:`'lr'` (left-to-right),
            :code:`'tb'` (top-to-bottom), :code:`'rl'` (right-to-left)
            (default is :code:`'bt'`).
        bgcolor : multiple types
            The background colour (default is :code:`'white'`).
        fgcolor : multiple types
            The foreground colour, given by default to nodes
            and labels (default is :code:`'black'`).
        labels : :class:`bool`
            Whether to display node labels (default is :code:`True`).
        inputcolor : multiple types
            The colour of edges corresponding to input faces
            (default is :code:`'magenta'`).
        outputcolor : multiple types
            The colour of edges corresponding to output faces
            (default is :code:`'blue'`).
        xscale : :class:`float`
            (TikZ only) Scale factor to apply to x axis in output
            (default is based on the dimension and maximal number of
            elements in one dimension).
        yscale : :class:`float`
            (TikZ only) Scale factor to apply to y axis in output
            (default is based on the dimension and maximal number of
            elements in one dimension).
        """

        # Parameters
        tikz = params.get('tikz', DEFAULT['tikz'])
        show = params.get('show', DEFAULT['show'])
        path = params.get('path', None)

        orientation = params.get(
                'orientation', DEFAULT['orientation'])

        xscale = params.get('xscale', None)
        yscale = params.get('yscale', None)
        SCALE = (
                max(
                    (len(self.nodes[n]) for n in range(self.nodes.dim+1)),
                    default=0),
                2*self.nodes.dim
                )
        if orientation in ('bt', 'tb'):
            xscale = SCALE[0] if xscale is None else xscale
            yscale = SCALE[1] if yscale is None else yscale
        if orientation in ('lr', 'rl'):
            xscale = SCALE[1] if xscale is None else xscale
            yscale = SCALE[0] if yscale is None else yscale

        bgcolor = params.get(
                'bgcolor', DEFAULT['bgcolor'])
        fgcolor = params.get(
                'fgcolor', DEFAULT['fgcolor'])

        labels = params.get(
                'labels', DEFAULT['labels'])

        inputcolor = params.get(
                'inputcolor', DEFAULT['inputcolor'])
        outputcolor = params.get(
                'outputcolor', DEFAULT['outputcolor'])

        coord = self.place_nodes()

        backendclass = drawing.TikZBackend if tikz else drawing.MatBackend
        backend = backendclass(
                bgcolor=bgcolor,
                fgcolor=fgcolor,
                orientation=orientation)

        for node in self.nodes:
            if labels:
                label = '{},{}'.format(
                        node.pos,
                        self.labels[node.dim][node.pos])
            else:
                label = node.pos
            backend.draw_label(
                label,
                coord[node],
                (0, 0),
                ha='center',
                va='center')
        for edge in self.diagram.edges:
            color = inputcolor if self.diagram.edges[edge]['sign'] == '-' \
                    else outputcolor
            backend.draw_arrow(
                coord[edge[0]],
                coord[edge[1]],
                color=color,
                shorten=0.8)

        backend.output(path=path, show=show, xscale=xscale, yscale=yscale)


def draw(*ogps, **params):
    """
    Given any number of oriented graded posets, or maps, or diagrams,
    generates their Hasse diagrams and draws them.

    This is the same as generating the Hasse diagram for each
    argument, and calling :meth:`Hasse.draw` with the given
    parameters on each one of them.

    Arguments
    ---------
    *ogps : :class:`ogposets.OgPoset | ogposets.OgMap | diagrams.Diagram`
        Any number of oriented graded posets or maps or diagrams.

    Keyword arguments
    -----------------
    **params
        Passed to :meth:`Hasse.draw`.
    """
    for ogp in ogps:
        Hasse(ogp).draw(**params)
