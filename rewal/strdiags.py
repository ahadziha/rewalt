"""
Implements string diagram visualisations.
"""

import networkx as nx

from rewal import utils
from rewal.diagrams import Diagram


class StrDiag:
    """
    Class for string diagrams.
    """
    def __init__(self, diagram, **kwargs):
        utils.typecheck(diagram, {'type': Diagram})

        base_graph = nx.DiGraph()
        if diagram.dim > 0:
            for x in diagram.shape[-1]:
                isdegenerate = x.dim == diagram.ambient[diagram[x]].dim
                base_graph.add_node(
                        x,
                        type='node',
                        degenerate=isdegenerate)
