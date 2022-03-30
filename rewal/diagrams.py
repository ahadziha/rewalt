"""
Implements diagrammatic sets and diagrams.
"""

from rewal import utils
from rewal.shapes import (Shape, ShapeMap)


class DiagSet:
    """
    Class for diagrammatic sets.
    """
    def __init__(self, **kwargs):
        """ Initialises to an empty diagrammatic set. """
        self._generators = dict()
        self._by_dim = dict()

        for key, value in kwargs:  # TODO: accepted kwargs
            setattr(self, key, value)

    def __eq__(self, other):
        return isinstance(other, DiagSet) and \
                self.generators == other.generators

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
        else:
            utils.typecheck(input, {'type': Diagram})

        if output is None:
            output = Diagram(self)
        else:
            utils.typecheck(output, {'type': Diagram})

        boundary = Shape.atom_cospan(
                input.shape, output.shape)
        shape = boundary.target
        mapping = [
                [None for _ in n_data]
                for n_data in shape.face_data
                ]

        for x in input.shape:
            y = boundary.fst[x]
            mapping[y.dim][y.pos] = input[x].name
        for x in output.shape:
            y = boundary.snd[x]
            if mapping[y.dim][y.pos] is None:
                mapping[y.dim][y.pos] = output[x].name
            else:
                if mapping[y.dim][y.pos] != output[x].name:
                    raise ValueError(utils.value_err(
                            output, 'boundary does not match '
                            'boundary of {}'.format(repr(input))))
        mapping[-1][0] = name

        new = Diagram._new(
                shape,
                self,
                mapping,
                name)

        # TODO: valid kwargs
        self._generators.update({
                name: {
                    'cell': new,
                    'faces': set(),
                    'cofaces': set(),
                    **kwargs
                    }
                })

        if new.dim in self._by_dim:
            self._by_dim[new.dim].add(name)
        else:
            self._by_dim[new.dim] = {name}

        if new.dim > 0:
            for x in mapping[-2]:
                self.generators[x]['cofaces'].add(name)
                self.generators[name]['faces'].add(x)

    def remove(self, name):
        """
        Removes a generator.
        """
        to_remove = [*self.generators[name]['cofaces']]
        for x in to_remove:
            self.remove(x)

        dim = self[name].dim
        self._by_dim[dim].remove(name)
        if len(self._by_dim[dim]) == 0:
            self._by_dim.pop(dim, None)

        for x in self.generators[name]['faces']:
            self.generators[x]['cofaces'].remove(name)

        self._generators.pop(name, None)


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
        self._name = ''

    def __eq__(self, other):
        return isinstance(other, Diagram) and \
                self.shape == other.shape and \
                self.ambient == other.ambient and \
                self.mapping == other.mapping

    def __getitem__(self, element):
        if element in self.shape:
            return self.ambient[self.mapping[element.dim][element.pos]]
        raise ValueError(utils.value_err(
            element, 'not an element of the shape'))

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

    @property
    def isdegenerate(self):
        if self.dim <= 0:
            return False
        for x in self.mapping[-1]:
            if self.ambient[x].dim == self.dim:
                return False
        return True

    def paste(self, other, dim=None):
        utils.typecheck(other, {
                'type': Diagram,
                'st': lambda x: x.ambient == self.ambient,
                'why': 'not the same ambient DiagSet'})
        if dim is None:
            dim = min(self.dim, other.dim) - 1

        paste_cospan = Shape.paste_cospan(
                self.shape, other.shape, dim)
        shape = paste_cospan.target
        mapping = [
                [None for _ in n_data]
                for n_data in shape.face_data
                ]

        for x in self.shape:
            y = paste_cospan.fst[x]
            mapping[y.dim][y.pos] = self[x].name
        for x in other.shape:
            y = paste_cospan.snd[x]
            if mapping[y.dim][y.pos] is None:
                mapping[y.dim][y.pos] = other[x].name
            else:
                if mapping[y.dim][y.pos] != other[x].name:
                    raise ValueError(utils.value_err(
                            other, 'input {}-boundary does not match '
                            'output {}-boundary of {}'.format(
                                str(dim), str(dim), repr(self))))

        name = '({} #{} {})'.format(
                str(self.name), str(dim), str(other.name))

        pasted = Diagram._new(
                shape,
                self.ambient,
                mapping,
                name)
        return pasted

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
                    self[x].name
                    for x in n_data
                ]
                for n_data in shapemap.mapping]
        return Diagram._new(
                shape,
                self.ambient,
                mapping,
                name)

    def boundary(self, sign, dim=None):
        """
        Boundaries of the diagram.
        """
        if dim is None:
            dim = self.dim - 1
        sign = utils.mksign(sign)
        name = '∂[{},{}]{}'.format(
                sign, str(dim), str(self.name))
        return self.pullback(self.shape.boundary_inclusion(
            sign, dim), name)

    def unit(self):
        """
        Unit on the diagram.
        """
        return self.pullback(self.shape.inflate(), self.name)

    def unitor_l(self, sign):
        """
        Left unitor or inverse unitor.
        """
        sign = utils.mksign(sign)
        unitor_map = self.shape.inflate(
                self.shape.boundary('+'))
        if sign == '-':
            unitor_map = unitor_map.dual(self.dim + 1)
        return self.pullback(unitor_map, self.name)

    def unitor_r(self, sign):
        """
        Right unitor or inverse unitor.
        """
        sign = utils.mksign(sign)
        unitor_map = self.shape.inflate(
                self.shape.boundary('-'))
        if sign == '+':
            unitor_map = unitor_map.dual(self.dim + 1)
        return self.pullback(unitor_map, self.name)

    def unitor(self, collapsed, invert=False):
        """
        Generic unitor.
        """
        unitor_map = self.shape.inflate(
                collapsed)
        if invert:
            unitor_map = unitor_map.dual(self.dim + 1)
        return self.pullback(unitor_map, self.name)

    # Internal methods
    @classmethod
    def _new(cls, shape, ambient, mapping, name=None):
        new = Diagram.__new__(cls)
        new._shape = shape
        new._ambient = ambient
        new._mapping = mapping
        new._name = name
        return new
