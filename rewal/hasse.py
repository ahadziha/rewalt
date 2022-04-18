"""
Implements oriented Hasse diagram visualisation.
"""

import networkx as nx

import rewal
from rewal import utils
from rewal.drawing import MatBackend


class Hasse:
    """
    Class for oriented Hasse diagrams.
    """
    def __init__(self, ogp):
        if isinstance(ogp, rewal.ogposets.OgPoset):
            self._labels = ogp.id().mapping
        elif isinstance(ogp, rewal.diagrams.Diagram):
            self._labels = ogp.mapping
            ogp = ogp.shape
        elif isinstance(ogp, rewal.ogposets.OgMap):
            self._labels = ogp.mapping
            ogp = ogp.source
        else:
            raise TypeError(utils.type_err(
                rewal.ogposets.OgPoset, ogp))

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
        ystep = 1/(dim+2)
        xstep = [1/(len(self.nodes[n])+1) for n in range(dim+1)]

        coordinates = dict()
        for x in self.nodes:
            coordinates[x] = (
                    (x.pos + 1)*xstep[x.dim],
                    (x.dim + 1)*ystep
                    )
        return coordinates

    def draw(self, **params):
        coord = self.place_nodes()

        backend = MatBackend()

        for node in self.nodes:
            backend.draw_label(
                node.pos,
                coord[node],
                (0, 0),
                ha='center',
                va='center')
        for edge in self.diagram.edges:
            color = 'magenta' if self.diagram.edges[edge]['sign'] == '-' \
                    else 'blue'
            backend.draw_arrow(
                coord[edge[0]],
                coord[edge[1]],
                color=color)

        backend.output()
