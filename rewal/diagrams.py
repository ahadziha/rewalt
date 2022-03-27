"""
Implements diagrammatic sets and diagrams.
"""

from rewal import utils
from rewal.shapes import Shape


class DiagSet:
    """
    Class for diagrammatic sets.
    """
    def __init__(self):
        """ Initialises to an empty diagrammatic set. """
        self._generators = dict()
        self._by_dim = dict()
        self._saved = dict()

    def __getitem__(self, key):
        if key in self.generators:
            return self._generators[key]['cell']
        raise KeyError(str(key))

    @property
    def generators(self):
        return self._generators.keys()

    def new_point(self, name, **kwargs):
        """ Adds a 0-cell. """
        if name in self.generators:
            raise ValueError(utils.value_err(
                name, 'already in use'))

        new = Diagram._new(
                Shape.point(),
                self,
                [name],
                name)

        self._generators.update({
                name: {
                    'cell': new,
                    **kwargs
                    }
                })
        if 0 in self._by_dim:
            self._by_dim[0].add(name)
        else:
            self._by_dim[0] = {name}


class Diagram:
    """
    Class for diagrams in diagrammatic sets.
    """
    def __init__(self, ambient):
        """ Initialises to the empty diagram. """
        utils.typecheck(ambient, {'type': DiagSet})

        self._shape = Shape.empty()
        self._ambient = ambient
        self._mapping = []
        self._name = 'empty'

    @property
    def name(self):
        return self._name

    def save(self, name=None):
        """
        Adds the current diagram to the list of saved diagrams.
        """
        if name is None:
            name = self.name

        self.ambient._saved.update(
                {
                    name: self
                })

    # Internal methods
    @classmethod
    def _new(cls, shape, ambient, mapping, name):
        new = Diagram.__new__(cls)
        new._shape = shape
        new._ambient = ambient
        new._mapping = mapping
        new._name = name
