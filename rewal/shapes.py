"""
Implements oriented graded posets and molecules.
"""

from rewal import utils


class El:
    """
    Defines an element of an oriented graded poset.
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
    Defines an oriented graded poset.
    """

    @staticmethod
    def __wfcheck(face_data):
        """
        Checks that face data is well formed.
        """

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
            # Check that faces are within bounds
            k = max([i for x in n_data for a in x for i in x[a]])
            if n == 0 or k >= sizes[n-1]:
                raise ValueError(utils.value_err(k, 'out of bounds'))

            # Check that input/output are inhabited and disjoint
            for x in n_data:
                if not x['-'].isdisjoint(x['+']):
                    raise ValueError(
                            'Input and output faces must be disjoint.')
                if n > 0 and x['-'] == x['+'] == set():
                    raise ValueError(
                            'The element must have at least one face.')

    def __init__(self, face_data, coface_data,
                 wfcheck=True, matchcheck=True):
        if wfcheck:
            self.__wfcheck(face_data)

        if matchcheck:
            pass
            # raise ValueError('Face and coface data do not match.')

        self._face_data = face_data
        self._coface_data = coface_data

    @property
    def face_data(self):
        return self._face_data

    @property
    def coface_data(self):
        return self._coface_data

    # TODO: __getitem__

    @property
    def size(self):
        return [len(_) for _ in self.face_data]

    def __eq__(self, other):
        return self.face_data == other.face_data

    @classmethod
    def from_face_data(cls, face_data, wfcheck=True):
        coface_data = [
                [
                    {'-': set(), '+': set()}
                    for _ in face_data[n]
                ] 
                for n in range(len(face_data))]
        for n, sn_elements in enumerate(face_data[1:]):
            pass
