"""
Implements diagrammatic sets and diagrams.
"""

import rewal
from rewal import utils
from rewal.shapes import (Shape, ShapeMap)


class DiagSet:
    """
    Class for diagrammatic sets.
    """
    _PRIVATE = (
            'shape',
            'mapping',
            'faces',
            'cofaces',
            'inverse')

    def __init__(self):
        """ Initialises to an empty diagrammatic set. """
        self._generators = dict()
        self._by_dim = dict()

    def __str__(self):
        return '{} with {} generators'.format(
                type(self).__name__, str(len(self)))

    def __eq__(self, other):
        return isinstance(other, DiagSet) and \
                self.generators == other.generators

    def __getitem__(self, key):
        if key in self.generators:
            return Diagram._new(
                    self.generators[key]['shape'],
                    self,
                    self.generators[key]['mapping'],
                    key)
        raise KeyError(str(key))

    def __len__(self):
        return len(self.generators)

    def __iter__(self):
        return iter(self.generators)

    @property
    def generators(self):
        return self._generators

    @property
    def by_dim(self):
        return self._by_dim

    @property
    def issimplicial(self):
        for x in self:
            if not isinstance(self[x], SimplexDiagram):
                return False
        return True

    @property
    def iscubical(self):
        for x in self:
            if not isinstance(self[x], CubeDiagram):
                return False
        return True

    def add(self, name, input=None, output=None, **kwargs):
        """ Adds a generator and returns it. """
        if name in self.generators:
            raise ValueError(utils.value_err(
                name, 'name already in use'))
        for key in self._PRIVATE:
            kwargs.pop(key, None)

        if input is None:
            input = Diagram(self)
        if output is None:
            output = Diagram(self)
        for x in (input, output):
            utils.typecheck(x, {
                'type': Diagram,
                'st': lambda x: x.ambient == self,
                'why': 'not a diagram in {}'.format(repr(self))})

        boundary = Shape.atom(
                input.shape, output.shape, cospan=True)
        shape = boundary.target
        mapping = [
                [None for _ in n_data]
                for n_data in shape.face_data
                ]

        for x in input.shape:
            y = boundary.fst[x]
            mapping[y.dim][y.pos] = input.mapping[x.dim][x.pos]
        for x in output.shape:
            y = boundary.snd[x]
            if mapping[y.dim][y.pos] is None:
                mapping[y.dim][y.pos] = output.mapping[x.dim][x.pos]
            else:
                if mapping[y.dim][y.pos] != output.mapping[x.dim][x.pos]:
                    raise ValueError(utils.value_err(
                            output, 'boundary does not match '
                            'boundary of {}'.format(repr(input))))
        mapping[-1][0] = name
        self._update_generators(name, shape, mapping, **kwargs)

        return self[name]

    def add_simplex(self, name, *faces, **kwargs):
        if len(faces) <= 1:
            return self.add(name, **kwargs)

        if name in self.generators:
            raise ValueError(utils.value_err(
                name, 'name already in use'))

        dim = len(faces) - 1
        for x in faces:
            utils.typecheck(x, {
                'type': SimplexDiagram,
                'st': lambda x: x.ambient == self and x.dim == dim-1,
                'why': 'expecting a {}-simplex in {}'.format(
                    str(dim - 1), repr(self))})

        shape = Shape.simplex(dim)
        mapping = [
                [None for _ in n_data]
                for n_data in shape.face_data
                ]

        for k, face in enumerate(faces):
            face_map = shape.simplex_face(k)
            for x in face:
                y = face_map[x]
                if mapping[y.dim][y.pos] is None:
                    mapping[y.dim][y.pos] = face.mapping[x.dim][x.pos]
                else:
                    if mapping[y.dim][y.pos] != face.mapping[x.dim][x.pos]:
                        raise ValueError(utils.value_err(
                            face, 'boundary of face does not '
                            'match other faces'))

        mapping[-1][0] = name
        self._update_generators(name, shape, mapping, **kwargs)

        return self[name]

    def add_cube(self, name, *faces, **kwargs):
        if len(faces) % 2 == 1:
            raise ValueError(utils.value_err(
                faces, 'expecting an even number of faces'))
        if name in self.generators:
            raise ValueError(utils.value_err(
                name, 'name already in use'))

        dim = int(len(faces) / 2)
        for x in faces:
            utils.typecheck(x, {
                'type': CubeDiagram,
                'st': lambda x: x.ambient == self and x.dim == dim-1,
                'why': 'expecting a {}-cube in {}'.format(
                    str(dim - 1), repr(self))})

        shape = Shape.cube(dim)
        mapping = [
                [None for _ in n_data]
                for n_data in shape.face_data
                ]

        for n, face in enumerate(faces):
            k = int(n/2)
            sign = utils.mksign(n % 2)
            face_map = shape.cube_face(k, sign)
            for x in face:
                y = face_map[x]
                if mapping[y.dim][y.pos] is None:
                    mapping[y.dim][y.pos] = face.mapping[x.dim][x.pos]
                else:
                    if mapping[y.dim][y.pos] != face.mapping[x.dim][x.pos]:
                        raise ValueError(utils.value_err(
                            face, 'boundary of face does not '
                            'match other faces'))

        mapping[-1][0] = name
        self._update_generators(name, shape, mapping, **kwargs)

        return self[name]

    def invert(self, name,
               inversename=None,
               rightinvertorname=None,
               leftinvertorname=None,
               **kwargs):
        """
        Adds an inverse and 'invertors' for a generator.
        """
        generator = self[name]
        if generator.dim == 0:
            raise ValueError(utils.value_err(
                name, 'cannot invert 0-cell'))
        if 'inverse' in self.generators[name]:
            raise ValueError(utils.value_err(
                name, 'already inverted'))

        if inversename is None:
            inversename = '{}⁻¹'.format(str(name))
        if rightinvertorname is None:
            rightinvertorname = 'inv({}, {})'.format(
                    str(name), str(inversename))
        if leftinvertorname is None:
            leftinvertorname = 'inv({}, {})'.format(
                    str(inversename), str(name))

        self._generators[name].update(
                {'inverse': inversename})

        inverse = self.add(
                inversename,
                generator.output,
                generator.input,
                inverse=name,
                **kwargs)
        leftinvertor = self.add(
                rightinvertorname,
                generator.paste(inverse),
                generator.input.unit())
        rightinvertor = self.add(
                leftinvertorname,
                inverse.paste(generator),
                generator.output.unit())

        return inverse, leftinvertor, rightinvertor

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

        if 'inverse' in self.generators[name]:
            inverse = self.generators[name]['inverse']
            self.generators[inverse].pop('inverse', None)

        for x in self.generators[name]['faces']:
            self.generators[x]['cofaces'].remove(name)

        self.generators.pop(name, None)

    def update(self, name, **kwargs):
        """
        Updates the optional arguments of a generator.
        """
        for key in self._PRIVATE:
            if key in kwargs:
                raise AttributeError(key, 'private attribute')
        self.generators[name].update(kwargs)

    def copy(self):
        new = DiagSet()
        new._generators = self.generators.copy()
        new._by_dim = self.by_dim.copy()
        return new

    @staticmethod
    def yoneda(shape):
        utils.typecheck(shape, {'type': Shape})
        yoneda = DiagSet()
        for x in shape:
            atom = shape.atom_inclusion(x)
            yoneda._generators.update({
                x: {
                    'shape': atom.source,
                    'mapping': atom.mapping,
                    'faces': {y for y in shape.faces(x)},
                    'cofaces': {y for y in shape.cofaces(x)}
                    }
                })
        for n, n_size in enumerate(shape.size):
            yoneda._by_dim.update({
                n: {
                    rewal.ogposets.El(n, k) for k in range(n_size)
                    }
                })
        return yoneda

    # Internal methods
    def _update_generators(self, name, shape, mapping, **kwargs):
        self._generators.update({
                name: {
                    'shape': shape,
                    'mapping': mapping,
                    'faces': set(),
                    'cofaces': set(),
                    **kwargs
                    }
                })

        if shape.dim in self._by_dim:
            self._by_dim[shape.dim].add(name)
        else:
            self._by_dim[shape.dim] = {name}

        if shape.dim > 0:
            for x in mapping[-2]:
                self.generators[x]['cofaces'].add(name)
                self.generators[name]['faces'].add(x)


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

    def __str__(self):
        return str(self.name)

    def __eq__(self, other):
        return isinstance(other, Diagram) and \
                self.shape == other.shape and \
                self.ambient == other.ambient and \
                self.mapping == other.mapping

    def __len__(self):
        return len(self.shape)

    def __getitem__(self, element):
        if element in self.shape:
            return self.ambient[self.mapping[element.dim][element.pos]]
        raise ValueError(utils.value_err(
            element, 'not an element of the shape'))

    def __iter__(self):
        return iter(self.shape)

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

    @property
    def iscell(self):
        """
        Returns whether the diagram is a cell, that is, its shape
        is an atom.
        """
        return self.shape.isatom

    def rename(self, name):
        """ Renames the diagram. """
        self._name = name

    # Methods for creating new diagrams
    def paste(self, other, dim=None,
              cospan=False):
        utils.typecheck(other, {
                'type': Diagram,
                'st': lambda x: x.ambient == self.ambient,
                'why': 'not the same ambient DiagSet'})
        if dim is None:
            dim = min(self.dim, other.dim) - 1

        paste_cospan = Shape.paste(
                self.shape, other.shape, dim, cospan=True)

        shape = paste_cospan.target
        mapping = Diagram._paste_fill_mapping(
                self, other, paste_cospan)
        name = '({} #{} {})'.format(
                str(self.name), str(dim), str(other.name))

        pasted = Diagram._new(
                shape,
                self.ambient,
                mapping,
                name)
        if cospan:
            return pasted, paste_cospan
        return pasted

    def paste_along(self, other, fst, snd, name=None,
                    cospan=False):
        utils.typecheck(other, {
                'type': Diagram,
                'st': lambda x: x.ambient == self.ambient,
                'why': 'not the same ambient DiagSet'})
        assoc = ((fst, self), (snd, other))
        for x in assoc:
            utils.typecheck(x[0], {
                'type': ShapeMap,
                'st': lambda f: f.target == x[1].shape,
                'why': 'expecting a map with target {}'.format(
                    repr(x[1].shape))}
                )
        paste_cospan = Shape.paste_along(
                fst, snd, cospan=True)

        shape = paste_cospan.target
        mapping = Diagram._paste_fill_mapping(
                self, other, paste_cospan)
        if name is None:
            name = '({} ## {})'.format(
                    str(self.name), str(other.name))

        pasted = Diagram._new(
                shape,
                self.ambient,
                mapping,
                name)
        if cospan:
            return pasted, paste_cospan
        return pasted

    def to_output(self, pos, other,
                  cospan=False):
        name = '([{} -> {}]{})'.format(
                str(other.name), str(pos), str(self.name))
        return self.paste_along(
                other,
                self.shape.atom_inclusion(
                    rewal.ogposets.El(self.dim-1, pos)),
                other.shape.boundary('-'),
                name,
                cospan)

    def to_input(self, pos, other,
                 cospan=False):
        name = '({}[{} <- {}])'.format(
                str(self.name), str(pos), str(other.name))
        return other.paste_along(
                self,
                other.shape.boundary('+'),
                self.shape.atom_inclusion(
                    rewal.ogposets.El(self.dim-1, pos)),
                name,
                cospan)

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
                    self.mapping[x.dim][x.pos]
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
        return self.pullback(self.shape.boundary(
            sign, dim), name)

    @property
    def input(self):
        return self.boundary('-')

    @property
    def output(self):
        return self.boundary('+')

    def unit(self):
        """
        Unit on the diagram.
        """
        name = '1{}'.format(str(self.name))
        return self.pullback(self.shape.inflate(), name)

    def unitor_l(self, sign):
        """
        Left unitor or inverse unitor.
        """
        sign = utils.mksign(sign)
        unitor_map = self.shape.inflate(
                self.shape.all().boundary('+'))
        if sign == '-':
            unitor_map = unitor_map.dual(self.dim + 1)
        name = 'λ[{}]{}'.format(
                sign, str(self.name))

        return self.pullback(unitor_map, name)

    def unitor_r(self, sign):
        """
        Right unitor or inverse unitor.
        """
        sign = utils.mksign(sign)
        unitor_map = self.shape.inflate(
                self.shape.all().boundary('-'))
        if sign == '+':
            unitor_map = unitor_map.dual(self.dim + 1)
        name = 'ρ[{}]{}'.format(
                sign, str(self.name))

        return self.pullback(unitor_map, name)

    def unitor(self, collapsed, invert=False):
        """
        Generic unitor.
        """
        unitor_map = self.shape.inflate(
                collapsed)
        if invert:
            unitor_map = unitor_map.dual(self.dim + 1)

        name = 'unitor on {}'.format(str(self.name))
        return self.pullback(unitor_map, name)

    # Alternative constructors
    @staticmethod
    def yoneda(shapemap, name=None):
        return Diagram._new(
                shapemap.source,
                DiagSet.yoneda(shapemap.target),
                shapemap.mapping,
                name)

    @staticmethod
    def with_layers(fst, *layers):
        """ Alternative constructor for LayeredDiagram. """
        return LayeredDiagram(fst, *layers)

    # Internal methods
    @staticmethod
    def _new(shape, ambient, mapping, name=None):
        def diagramclass():
            if isinstance(shape, rewal.shapes.Point):
                return PointDiagram
            if isinstance(shape, rewal.shapes.Arrow):
                return ArrowDiagram
            if isinstance(shape, rewal.shapes.Cube):
                return CubeDiagram
            if isinstance(shape, rewal.shapes.Simplex):
                return SimplexDiagram
            return Diagram

        new = Diagram.__new__(diagramclass())
        new._shape = shape
        new._ambient = ambient
        new._mapping = mapping
        new._name = name
        return new

    @staticmethod
    def _paste_fill_mapping(fst, snd, paste_cospan):
        shape = paste_cospan.target
        mapping = [
            [None for _ in n_data]
            for n_data in shape.face_data
            ]

        for x in fst.shape:
            y = paste_cospan.fst[x]
            mapping[y.dim][y.pos] = fst.mapping[x.dim][x.pos]
        for x in snd.shape:
            y = paste_cospan.snd[x]
            if mapping[y.dim][y.pos] is None:
                mapping[y.dim][y.pos] = snd.mapping[x.dim][x.pos]
            else:
                if mapping[y.dim][y.pos] != snd.mapping[x.dim][x.pos]:
                    raise ValueError(utils.value_err(
                            snd, 'boundary does not match '
                            'boundary of {}'.format(repr(fst))))
        return mapping


class SimplexDiagram(Diagram):
    def simplex_face(self, k):
        face_map = self.shape.simplex_face(k)
        name = 'd[{}]{}'.format(
                str(k), str(self.name))
        return self.pullback(face_map, name)

    def simplex_degeneracy(self, k):
        degeneracy_map = self.shape.simplex_degeneracy(k)
        name = 's[{}]{}'.format(
                str(k), str(self.name))
        return self.pullback(degeneracy_map, name)


class CubeDiagram(Diagram):
    def cube_face(self, k, sign):
        sign = utils.mksign(sign)
        face_map = self.shape.cube_face(k, sign)
        name = 'δ[{},{}]{}'.format(
                str(k), sign, str(self.name))
        return self.pullback(face_map, name)

    def cube_degeneracy(self, k):
        degeneracy_map = self.shape.cube_degeneracy(k)
        name = 'σ[{}]{}'.format(
                str(k), str(self.name))
        return self.pullback(degeneracy_map, name)

    def cube_connection(self, k, sign):
        sign = utils.mksign(sign)
        connection_map = self.shape.cube_connection(k, sign)
        name = 'γ[{},{}]{}'.format(
                str(k), sign, str(self.name))
        return self.pullback(connection_map, name)


class ArrowDiagram(SimplexDiagram, CubeDiagram):
    pass


class PointDiagram(SimplexDiagram, CubeDiagram):
    def degeneracy(self, shape):
        utils.typecheck(shape, {'type': rewal.shapes.Shape})
        return self.pullback(
                shape.terminal(), self.name)


class LayeredDiagram(Diagram):
    """
    A diagram that remembers a layering in the top dimension.
    """
    def __init__(self, fst, *others):
        utils.typecheck(fst, {'type': Diagram})
        dim = fst.dim

        diagram = fst
        layers = [fst.shape.id()]
        for x in others:
            utils.typecheck(x, {
                'type': Diagram,
                'st': lambda x: x.dim == dim,
                'why': 'expecting diagram of dimension {}'.format(
                    str(dim))})
            diagram, cospan = diagram.paste(x, cospan=True)
            layers = [
                    *[f.then(cospan.fst) for f in layers],
                    cospan.snd]

        self._shape = diagram.shape
        self._ambient = diagram.ambient
        self._mapping = diagram.mapping
        self._name = diagram.name

        self._layers = layers

    @property
    def layers(self):
        return [
                self.pullback(f, 'layer {} of {}'.format(
                    str(n), self.name))
                for n, f in enumerate(self._layers)]

    @property
    def rewrite_steps(self):
        rewrite_steps = [
                *[layer.input for layer in self.layers],
                self.layers[-1].output
                ]
        for n, step in enumerate(rewrite_steps):
            step.rename('step {} of {}'.format(
                    str(n), self.name))
        return rewrite_steps

    def concatenate(self, *others):
        for x in others:
            utils.typecheck(x, {'type': LayeredDiagram})

        if len(others) == 0:
            return self
        if len(others) > 1:
            return self.concatenate(
                    others[0]).concatenate(others[1:])
        other = others[0]
        diagram, cospan = self.paste(other, cospan=True)

        concat = Diagram.__new__(LayeredDiagram)
        concat._shape = diagram.shape
        concat._mapping = diagram.mapping
        concat._ambient = diagram.ambient
        concat._name = diagram.name

        concat._layers = [
            *[f.then(cospan.fst) for f in self._layers],
            *[f.then(cospan.snd) for f in other._layers]]
        return concat
