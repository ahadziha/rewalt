"""
Implements oriented graded posets and molecules.
"""

import numpy as np

from rewal import utils


class El:
    """
    Class for elements of an oriented graded poset.
    """

    def __init__(self, dim, pos):
        for x in (dim, pos):
            utils.typecheck(x, {
                'type': int,
                'st': lambda x: x >= 0,
                'why': 'expecting non-negative integer'
                })
        self._dim = dim
        self._pos = pos

    def __repr__(self):
        return "El({}, {})".format(repr(self.dim), repr(self.pos))

    def __str__(self):
        return repr(self)

    def __eq__(self, other):
        return isinstance(other, El) and \
                self.dim == other.dim and self.pos == other.pos

    def __hash__(self):
        return hash(repr(self))

    @property
    def dim(self):
        """
        The dimension of an element is immutable.
        """
        return self._dim

    @property
    def pos(self):
        """
        The position of an element is immutable.
        """
        return self._pos

    def shift(self, k):
        utils.typecheck(k, {
            'type': int,
            'st': lambda x: self.pos + x >= 0,
            'why': 'shifted position must be non-negative'
            })
        return El(self.dim, self.pos + k)


class OgPoset:
    """
    Class for oriented graded posets.
    """

    def __init__(self, face_data, coface_data,
                 wfcheck=True, matchcheck=True):
        if wfcheck:
            self._wfcheck(face_data)

        if matchcheck:
            if not coface_data == self._coface_from_face(face_data):
                raise ValueError("Face and coface data do not match.")

        self._face_data = face_data
        self._coface_data = coface_data

    def __str__(self):
        return "OgPoset with {} elements".format(str(self.size))

    def __getitem__(self, key):
        """ Returns a GrSubset. """
        return self.all_elements[key]

    def __iter__(self):
        return iter(self.all_elements)

    def __eq__(self, other):
        return isinstance(other, OgPoset) and \
                self.face_data == other.face_data and \
                self.coface_data == other.coface_data

    @property
    def face_data(self):
        """ Face data are read-only. """
        return self._face_data

    @property
    def coface_data(self):
        """ Coface data are read-only. """
        return self._coface_data

    @property
    def size(self):
        """ Returns the number of elements in each dimension as a list. """
        return [len(_) for _ in self.face_data]

    @property
    def dim(self):
        """ Returns the dimension of the oriented graded poset. """
        return len(self.face_data) - 1

    @property
    def chain(self):
        """ Returns chain complex representation. """
        size, dim = self.size, self.dim
        chain = [
                np.zeros((size[i], size[i+1]), dtype=int) for i in range(dim)
                ]
        for n, n_data in enumerate(self.coface_data):
            for i, x in enumerate(n_data):
                for j in x['-']:
                    chain[n][i][j] = -1
                for j in x['+']:
                    chain[n][i][j] = 1

        return chain

    @property
    def all_elements(self):
        """ Returns the GrSubset of all elements. """
        return GrSubset(
                GrSet(*[El(n, k) for n in range(len(self.size))
                        for k in range(self.size[n])]),
                self, wfcheck=False, mkfresh=False)

    def faces(self, element, sign=None):
        """
        Returns the graded set of faces of an element with given sign.
        """
        if sign is None:
            return self.faces(element, '-').union(
                self.faces(element, '+'))
        sign = utils.make_sign(sign)
        utils.typecheck(element, {'type': El})
        n = element.dim
        k = element.pos
        if n <= self.dim and k <= self.size[n]:
            return GrSet(*[El(n-1, i)
                           for i in self.face_data[n][k][sign]])

    def cofaces(self, element, sign=None):
        """
        Returns the graded set of cofaces of an element with given sign.
        """
        if sign is None:
            return self.cofaces(element, '-').union(
                self.cofaces(element, '+'))
        sign = utils.make_sign(sign)
        utils.typecheck(element, {'type': El})
        n = element.dim
        k = element.pos
        if n <= self.dim and k <= self.size[n]:
            return GrSet(*[El(n+1, i)
                           for i in self.coface_data[n][k][sign]])

    # Class methods
    @classmethod
    def from_face_data(cls, face_data, wfcheck=True):
        if wfcheck:
            cls._wfcheck(face_data)
        coface_data = cls._coface_from_face(face_data)

        return cls(face_data, coface_data, wfcheck=False, matchcheck=False)

    # Internal methods
    @staticmethod
    def _wfcheck(face_data):
        """ Internal method checking that face data is well-formed. """

        utils.typecheck(face_data, {'type': list}, {
            'type': list,
            'st': lambda x: len(x) > 0,
            'why': 'expecting non-empty list'
            }, {
            'type': dict,
            'st': lambda x: x.keys() == {'-', '+'},
            'why': "expecting dict with keys '-', '+'"
            }, {'type': set}, {
            'type': int,
            'st': lambda x: x >= 0,
            'why': 'expecting non-negative int'})

        sizes = [len(_) for _ in face_data]
        for n, n_data in enumerate(face_data):
            # Check that faces are within bounds.
            k = max([i for x in n_data for sign in x for i in x[sign]],
                    default=-1)
            if (n == 0 and k >= 0) or k >= sizes[n-1]:
                raise ValueError(utils.value_err(k, 'out of bounds'))

            # Check that input/output are inhabited and disjoint.
            for x in n_data:
                if not x['-'].isdisjoint(x['+']):
                    raise ValueError(
                            'Input and output faces must be disjoint.')
                if n > 0 and x['-'] == x['+'] == set():
                    raise ValueError(
                            'The element must have at least one face.')

    @staticmethod
    def _coface_from_face(face_data):
        """
        Internal method constructing coface data from face data.
        Face data is presumed to be well-formed.
        """

        coface_data = [
                [
                    {'-': set(), '+': set()}
                    for _ in face_data[n]
                ]
                for n in range(len(face_data))]

        for n, sn_data in enumerate(face_data[1:]):
            for k, x in enumerate(sn_data):
                for a in ('-', '+'):
                    for i in x[a]:
                        coface_data[n][i][a].add(k)

        return coface_data


class GrSet:
    """
    Class for graded sets of elements of an oriented graded poset.
    """

    def __init__(self, *elements):
        self._elements = {}
        for x in elements:
            self.add(x)

    def __repr__(self):
        return "GrSet{}".format(repr(self.as_list))

    def __str__(self):
        return repr(self)

    def __contains__(self, item):
        if isinstance(item, El):
            if item.dim in self.grades:
                if item.pos in self._elements[item.dim]:
                    return True
        return False

    def __len__(self):
        """ Returns the total number of elements. """
        return len(self.as_list)

    def __iter__(self):
        """ The iterator is the iterator of the element list. """
        return iter(self.as_list)

    def __getitem__(self, key):
        if isinstance(key, int) and key >= -2:
            if key in self.grades:
                return GrSet(*[El(key, k) for k in self._elements[key]])
            return GrSet()
        if isinstance(key, slice):
            stop = key.stop if key.stop is not None else self.dim + 1
            indices = list(range(stop)[key])
            return GrSet().union(*[self[n] for n in indices])
        raise KeyError(str(key))

    def __eq__(self, other):
        return isinstance(other, GrSet) and \
                self._elements == other._elements

    @property
    def grades(self):
        """
        Returns the list of dimensions in which the graded set is not empty.
        """
        return sorted(self._elements.keys())

    @property
    def dim(self):
        """
        Returns the maximal dimension in which the graded set is not empty,
        or -1 if it is empty.
        """
        return max(self.grades, default=-1)

    @property
    def as_set(self):
        """ Returns a non-graded set containing the same elements. """
        return {El(n, k) for n in self._elements for k in self._elements[n]}

    @property
    def as_list(self):
        """ Returns the list of elements in increasing order. """
        return [El(n, k) for n in sorted(self._elements)
                for k in sorted(self._elements[n])]

    def add(self, element):
        """ Adds an element. """
        utils.typecheck(element, {'type': El})

        if element.dim in self.grades:
            self._elements[element.dim].add(element.pos)
        else:
            self._elements[element.dim] = {element.pos}

    def remove(self, element):
        """ Removes an element. """
        if element not in self:
            raise ValueError('{} not in graded set.'.format(repr(element)))

        self._elements[element.dim].remove(element.pos)
        if self._elements[element.dim] == set():
            del self._elements[element.dim]

    def union(self, *others):
        """
        Returns the union of the graded set with other graded sets.
        """
        for x in others:
            utils.typecheck(x, {'type': GrSet})

        union_as_set = self.as_set.union(*[x.as_set for x in others])
        return GrSet(*union_as_set)

    def intersection(self, *others):
        """
        Returns the intersection of the graded set with other graded sets.
        """
        for x in others:
            utils.typecheck(x, {'type': GrSet})

        intersection_as_set = self.as_set.intersection(
                *[x.as_set for x in others])
        return GrSet(*intersection_as_set)

    def issubset(self, other):
        """ Returns whether the graded set is a subset of another. """
        utils.typecheck(other, {'type': GrSet})
        return self.as_set.issubset(other.as_set)

    def isdisjoint(self, other):
        """ Returns whether the graded set is disjoint from another. """
        utils.typecheck(other, {'type': GrSet})
        return self.as_set.isdisjoint(other.as_set)

    def copy(self):
        """ Returns a copy of the graded set. """
        return GrSet(*self)


class GrSubset:
    """
    Class for pairs of a GrSet and an "ambient" OgPoset, where the GrSet
    is seen as a subset of the ambient.
    """

    def __init__(self, graded_set, ambient,
                 wfcheck=True, mkfresh=True):
        if wfcheck:
            self._wfcheck(graded_set, ambient)

        self._graded_set = graded_set.copy() if mkfresh else graded_set
        self._ambient = ambient

    def __eq__(self, other):
        return isinstance(other, GrSubset) and \
                self.proj == other.proj and \
                self.ambient == other.ambient

    def __str__(self):
        return 'GrSubset with {} elements in {}'.format(
                str(len(self.proj)), str(self.ambient))

    def __contains__(self, item):
        return item in self.proj

    def __len__(self):
        return len(self.proj)

    def __getitem__(self, key):
        return GrSubset(self.proj[key], self.ambient,
                        wfcheck=False, mkfresh=False)

    def __iter__(self):
        return iter(self.proj)

    @property
    def proj(self):
        """ Returns the underlying graded set. """
        return self._graded_set

    @property
    def ambient(self):
        """ The ambient OgPoset is read-only. """
        return self._ambient

    def maximal(self, close_first=True):
        """
        Returns the subset of elements that are not below any other
        elements in the graded set.
        """
        maximal = GrSet()
        closed_self = self.closure() if close_first else self

        for x in closed_self:
            if self.ambient.cofaces(x).isdisjoint(
                    closed_self.proj[x.dim + 1]):
                maximal.add(x)
        return GrSubset(maximal, self.ambient,
                        wfcheck=False, mkfresh=False)

    def closure(self):
        """
        Returns the closure of the subset in the oriented graded poset.
        """
        closure = self.proj.copy()

        for n in range(self.proj.dim, 0, -1):
            for element in closure[n]:
                for face in self.ambient.faces(element):
                    closure.add(face)

        return GrSubset(closure, self.ambient,
                        wfcheck=False, mkfresh=False)

    @property
    def isclosed(self):
        """ Returns whether the subset is closed. """
        return self.closure() == self

    def boundary(self, sign, dim=None,
                 close_first=True):
        """
        Returns the input or output dim-boundary of (the closure of)
        the graded set.
        """
        sign = utils.make_sign(sign)
        if dim is None:
            dim = self.proj.dim - 1
        closed_self = self.closure() if close_first else self

        # Add lower-dim maximal elements
        boundary = self.maximal().proj[:dim]

        for x in closed_self[dim]:
            if self.ambient.cofaces(x, utils.flip(sign)).isdisjoint(
                    closed_self.proj[x.dim + 1]):
                boundary.add(x)

        return GrSubset(boundary, self.ambient,
                        wfcheck=False, mkfresh=False).closure()

    # Internal methods
    @staticmethod
    def _wfcheck(graded_set, ambient):
        utils.typecheck(graded_set, {'type': GrSet})
        utils.typecheck(ambient, {'type': OgPoset})

        if graded_set.dim > ambient.dim:
            raise ValueError('Not a valid graded subset.')
        for n in graded_set.grades:
            if max([x.pos for x in graded_set[n]]) > ambient.size[n]:
                raise ValueError('Not a valid graded subset.')


class OgMap:
    """
    Class for (partial) maps of oriented graded posets.
    """

    def __init__(self, source, target, mapping=None,
                 wfcheck=True):
        # TODO: whole thing needs rewriting.
        if wfcheck:
            self._wfcheck(source, target, mapping)

        self._source = source
        self._target = target
        self._mapping = mapping if mapping is not None \
            else [[None for _ in range(source.size[n])]
                  for n in range(source.dim + 1)]

    @property
    def source(self):
        return self._source

    @property
    def target(self):
        return self._target

    @property
    def mapping(self):
        return self._mapping

    def __eq__(self, other):
        return isinstance(other, OgMap) and \
                self.source == other.source and \
                self.target == other.target and \
                self.mapping == other.mapping

    # Internal methods.
    @staticmethod
    def _wfcheck(source, target, mapping):
        for x in (source, target):
            utils.typecheck(x, {'type': OgPoset})

        # TODO: all the actual checks...
        # Have this call “single extension” checks.
