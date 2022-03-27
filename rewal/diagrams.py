"""
Implements diagrammatic sets and diagrams.
"""

from rewal import utils
from rewal.shapes import (Shape, ShapeMap)


class DiagSet:
    """
    Class for diagrammatic sets.
    """
    def __init__(self):
        """ Initialises to an empty diagrammatic set. """
        self._generators = dict()
        self._by_dim = dict()
        self._saved = dict()
        self._coface_data = dict()

    def __getitem__(self, key):
        if key in self.generators:
            return self._generators[key]['cell']
        raise KeyError(str(key))

    @property
    def generators(self):
        return self._generators

    @property
    def by_dim(self):
        return self._by_dim

    def add(self, name, input=None, output=None, **kwargs):
        """ Adds a generator. """
        if name in self.generators:
            raise ValueError(utils.value_err(
                name, 'already in use'))
        if input is None:
            input = Diagram(self)
        if output is None:
            output = Diagram(self)
        for x in (input, output):
            utils.typecheck(x, {'type': Diagram})

        shape = Shape.atom(
                input.shape, output.shape)

        boundary_in = shape.boundary_inclusion('-')
        boundary_out = shape.boundary_inclusion('+')
        mapping = [
                [None for _ in n_data]
                for n_data in shape.face_data
                ]
        for x in input.shape:
            y = boundary_in[x]
            mapping[y.dim][y.pos] = input[x]
        for x in output.shape:
            y = boundary_out[x]
            if mapping[y.dim][y.pos] is None:
                mapping[y.dim][y.pos] = output[x]
            else:
                if mapping[y.dim][y.pos] != output[x]:
                    raise ValueError(utils.value_err(
                            output, 'boundary does not match '
                            'boundary of {}'.format(repr(input))))
        mapping[-1][0] = name

        new = Diagram._new(
                shape,
                self,
                mapping,
                name)

        self._generators.update({
                name: {
                    'cell': new,
                    **kwargs
                    }
                })

        if new.dim in self._by_dim:
            self._by_dim[new.dim].add(name)
        else:
            self._by_dim[new.dim] = {name}

        self._coface_data.update(
                {name: set()}
                )
        if new.dim > 0:
            for x in mapping[-2]:
                self._coface_data[x].add(name)


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

    def __eq__(self, other):
        return isinstance(other, Diagram) and \
                self.shape == other.shape and \
                self.ambient == other.ambient and \
                self.mapping == other.mapping

    def __getitem__(self, element):
        if element in self.shape:
            return self.mapping[element.dim][element.pos]
        raise ValueError(utils.value_err(
            element, 'not in shape'))

    @property
    def name(self):
        return self._name

    @property
    def shape(self):
        return self._shape

    @property
    def ambient(self):
        return self._ambient

    @property
    def mapping(self):
        return self._mapping

    @property
    def dim(self):
        return self.shape.dim

    def paste(self, other, dim=None, name=None):
        utils.typecheck(other, {
                'type': Diagram,
                'st': lambda x: x.ambient == self.ambient,
                'why': 'not the same ambient DiagSet'})
        if dim is None:
            dim = min(self.dim, other.dim) - 1

        paste_cospan = Shape._paste_cospan(
                self.shape, other.shape, dim)
        shape = paste_cospan.target
        mapping = [
                [None for _ in n_data]
                for n_data in shape.face_data
                ]

        for x in self.shape:
            y = paste_cospan.fst[x]
            mapping[y.dim][y.pos] = self[x]
        for x in other.shape:
            y = paste_cospan.snd[x]
            if mapping[y.dim][y.pos] is None:
                mapping[y.dim][y.pos] = other[x]
            else:
                if mapping[y.dim][y.pos] != other[x]:
                    raise ValueError(utils.value_err(
                            other, 'input {}-boundary does not match '
                            'output {}-boundary of {}'.format(
                                str(dim), str(dim), repr(self))))

        if name is None:
            name = '({} #{} {})'.format(
                    str(self.name), str(dim), str(other.name))

        pasted = Diagram._new(
                shape,
                self.ambient,
                mapping,
                name)
        return pasted

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

    def pullback(self, shapemap, name=None):
        """
        Pullback of the diagram by a ShapeMap.
        """
        utils.typecheck(shapemap, {
            'type': ShapeMap,
            'st': lambda x: x.target == self.shape,
            'why': 'target does not match diagram shape'})

        shape = shapemap.source
        mapping = [
                [
                    self[x]
                    for x in n_data
                ]
                for n_data in shapemap.mapping]
        return Diagram._new(
                shape,
                self.ambient,
                mapping,
                name)

    def boundary(self, sign, dim=None):
        return self.pullback(self.shape.boundary_inclusion(
            sign, dim))

    # Internal methods
    @classmethod
    def _new(cls, shape, ambient, mapping, name=None):
        new = Diagram.__new__(cls)
        new._shape = shape
        new._ambient = ambient
        new._mapping = mapping
        new._name = name
        return new
