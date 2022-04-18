"""
Implements oriented Hasse diagram visualisation.
"""

import networkx as nx

from rewal import utils, ogposets, diagrams, drawing


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
    Class for oriented Hasse diagrams.
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
        return self._nodes

    @property
    def diagram(self):
        return self._diagram

    @property
    def labels(self):
        return self._labels

    def place_nodes(self):
        dim = self.nodes.dim
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
        # Parameters
        tikz = params.get('tikz', DEFAULT['tikz'])
        show = params.get('show', DEFAULT['show'])
        path = params.get('path', None)

        xscale = params.get(
                'xscale',
                max(
                    (len(self.nodes[n]) for n in range(self.nodes.dim+1))
                    ))
        yscale = params.get(
                'yscale',
                2*self.nodes.dim)

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

        orientation = params.get(
                'orientation', DEFAULT['orientation'])

        coord = self.place_nodes()

        backendclass = drawing.TikZBackend if tikz else drawing.MatBackend
        backend = backendclass(
                bgcolor=bgcolor,
                fgcolor=fgcolor,
                orientation=orientation)

        for node in self.nodes:
            if labels:
                label = '{}({})'.format(
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
    for ogp in ogps:
        Hasse(ogp).draw(**params)
