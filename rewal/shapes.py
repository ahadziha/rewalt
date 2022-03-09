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

    def __repr__(self):
        return "El({}, {})".format(str(self.dim), str(self.pos))

    def __eq__(self, other):
        return self.dim == other.dim and self.pos == other.pos

    def __hash__(self):
        return hash(repr(self))

    def shift(self, k):
        utils.typecheck(k, {
            'type': int,
            'st': lambda x: self.pos + x >= 0,
            'why': 'out of bounds'
            })
        return El(self.dim, self.pos + k)


class OgPos:
    """
    Class for oriented graded posets.
    """

    def __init__(self, face_data, coface_data,
                 wfcheck=True, matchcheck=True):
        if wfcheck:
            self.__wfcheck(face_data)

        if matchcheck:
            if not coface_data == \
                    self.__coface_from_face(face_data):
                raise ValueError("Face and coface data do not match.")

        self._face_data = face_data
        self._coface_data = coface_data

    @property
    def face_data(self):
        return self._face_data

    @property
    def coface_data(self):
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

    # TODO: __repr__, __str__?

    def __getitem__(self, key):
        utils.typecheck(key, {
            'type': int,
            'st': lambda x: x in range(self.dim + 1),
            'why': 'out of bounds'
            })

        # TODO: may not be the most useful format
        return list(range(self.size[key]))

    def __eq__(self, other):
        return self.face_data == other.face_data

    # Class methods
    @classmethod
    def from_face_data(cls, face_data, wfcheck=True):
        if wfcheck:
            cls.__wfcheck(face_data)
        coface_data = cls.__coface_from_face(face_data)

        return cls(face_data, coface_data, wfcheck=False, matchcheck=False)

    # Private methods
    @staticmethod
    def __wfcheck(face_data):
        """ Private method checking that face data is well-formed. """

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
            k = max([i for x in n_data for a in x for i in x[a]], default=-1)
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
    def __coface_from_face(face_data):
        """
        Private method constructing coface data from face data.
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
