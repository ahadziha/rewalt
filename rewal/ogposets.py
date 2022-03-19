"""
Implements oriented graded posets, their subsets, and their maps.
"""

import numpy as np

from rewal import utils


class El(tuple):
    """
    Class for elements of an oriented graded poset.
    """

    def __new__(self, dim, pos):
        for x in dim, pos:
            utils.typecheck(x, {
                'type': int,
                'st': lambda n: n >= 0,
                'why': 'expecting non-negative integer'
                })
        return tuple.__new__(El, (dim, pos))

    def __repr__(self):
        return "El({}, {})".format(repr(self.dim), repr(self.pos))

    def __str__(self):
        return repr(self)

    def __eq__(self, other):
        return type(self) is type(other) and \
                self.dim == other.dim and self.pos == other.pos

    def __hash__(self):
        return hash(repr(self))

    @property
    def dim(self):
        return self[0]

    @property
    def pos(self):
        return self[1]

    def shifted(self, k):
        utils.typecheck(k, {
            'type': int,
            'st': lambda n: self.pos + n >= 0,
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
            OgPoset.__wfcheck(face_data)

        if matchcheck:
            if not coface_data == OgPoset.__coface_from_face(face_data):
                raise ValueError("Face and coface data do not match.")

        self._face_data = face_data
        self._coface_data = coface_data

    def __str__(self):
        return "{} with {} elements".format(
                type(self).__name__, str(self.size))

    def __getitem__(self, key):
        return self.all()[key]

    def __contains__(self, item):
        if isinstance(item, El):
            if item.dim <= self.dim:
                if item.pos < self.size[item.dim]:
                    return True
        return False

    def __len__(self):
        return sum(self.size)

    def __iter__(self):
        return iter(self.all())

    def __eq__(self, other):
        return type(self) is type(other) and \
                self.face_data == other.face_data and \
                self.coface_data == other.coface_data

    def __add__(self, other):
        return OgPoset.disjoint_union(self, other)

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
        """
        Returns the number of elements in each dimension as a list.
        """
        return [len(_) for _ in self.face_data]

    @property
    def dim(self):
        """
        Returns the dimension of the oriented graded poset.
        """
        return len(self.face_data) - 1

    @property
    def as_chain(self):
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

    def all(self):
        """ Returns the Closed subset of all elements. """
        return Closed(
                GrSet(*[El(n, k) for n in range(len(self.size))
                        for k in range(self.size[n])]),
                self,
                wfcheck=False)

    def none(self):
        """ Returns the empty Closed subset. """
        return Closed(GrSet(), self,
                      wfcheck=False)

    def maximal(self):
        """ Returns the GrSubset of maximal elements. """
        return self.all().maximal()

    def faces(self, element, sign=None):
        """
        Returns the faces of an element as a graded set.
        """
        if element.dim == 0:
            return GrSet()
        if sign is None:
            return self.faces(element, '-').union(
                self.faces(element, '+'))
        sign = utils.mksign(sign)
        utils.typecheck(element, {
            'type': El,
            'st': lambda x: x.dim <= self.dim and x.pos <= self.size[x.dim],
            'why': 'out of bounds'})
        return GrSet(
                *[El(element.dim - 1, i)
                  for i in self.face_data[element.dim][element.pos][sign]]
                )

    def cofaces(self, element, sign=None):
        """
        Returns the cofaces of an element as a graded set.
        """
        if element.dim == self.dim:
            return GrSet()
        if sign is None:
            return self.cofaces(element, '-').union(
                self.cofaces(element, '+'))
        sign = utils.mksign(sign)
        utils.typecheck(element, {
            'type': El,
            'st': lambda x: x.dim <= self.dim and x.pos <= self.size[x.dim],
            'why': 'out of bounds'})
        return GrSet(
                *[El(element.dim + 1, i)
                  for i in self.coface_data[element.dim][element.pos][sign]]
                )

    def id(self):
        """ Returns the identity map on the OgPoset. """
        mapping = [
                [El(n, k) for k in range(self.size[n])]
                for n in range(self.dim + 1)
                ]
        return OgMap(self, self, mapping,
                     wfcheck=False)

    def image(self, ogmap):
        """ Returns the image of the whole OgPoset through an OgMap. """
        return self.all().image(ogmap)

    def boundary_inclusion(self, sign=None, dim=None):
        """ Returns the inclusion of the n-boundary into the OgPoset. """
        if isinstance(dim, int) and dim >= self.dim:
            return self.id()
        return self.all().boundary(sign, dim).as_map

    def boundary(self, sign=None, dim=None):
        """ Returns the n-boundary as another OgPoset. """
        return self.boundary_inclusion(sign, dim).source

    @classmethod
    def from_face_data(cls, face_data,
                       wfcheck=True):
        if wfcheck:
            cls.__wfcheck(face_data)
        coface_data = cls.__coface_from_face(face_data)

        return cls(face_data, coface_data,
                   wfcheck=False, matchcheck=False)

    @staticmethod
    def empty():
        """ The initial oriented graded poset. """
        return OgPoset([], [],
                       wfcheck=False, matchcheck=False)

    @staticmethod
    def point():
        """ The terminal oriented graded poset. """
        return OgPoset(
            [[{'-': set(), '+': set()}]],
            [[{'-': set(), '+': set()}]],
            wfcheck=False, matchcheck=False)

    @staticmethod
    def coproduct(fst, snd):
        """ Returns the coproduct cospan of two OgPosets. """
        for x in fst, snd:
            utils.typecheck(x, {'type': OgPoset})

        # Need to ensure all the data have the same length.
        offset1 = max(snd.dim - fst.dim, 0)
        offset2 = max(fst.dim - snd.dim, 0)

        shift = fst.size
        shift = shift + [0 for _ in range(offset1)]

        face_data_fst = [
                [
                    {sign: faces for sign, faces in x.items()}
                    for x in n_data
                ]
                for n_data in fst.face_data] + [[] for _ in range(offset1)]
        face_data_snd = [
                [
                    {sign:
                        {k + shift[n-1] for k in faces}
                     for sign, faces in x.items()}
                    for x in n_data
                ]
                for n, n_data in enumerate(snd.face_data)] + [
                        [] for _ in range(offset2)]
        face_data = [x + y
                     for x, y in zip(face_data_fst, face_data_snd)]

        coface_data_fst = [
                [
                    {sign: cofaces for sign, cofaces in x.items()}
                    for x in n_data
                ]
                for n_data in fst.coface_data] + [
                        [] for _ in range(offset1)]
        coface_data_snd = [
                [
                   {sign:
                       {k + shift[n+1] for k in cofaces}
                    for sign, cofaces in x.items()}
                   for x in n_data
                ]
                for n, n_data in enumerate(snd.coface_data)] + [
                        [] for _ in range(offset2)]
        coface_data = [x + y
                       for x, y in zip(coface_data_fst, coface_data_snd)]

        disjoint_union = OgPoset(face_data, coface_data,
                                 wfcheck=False, matchcheck=False)

        mapping_fst = fst.id().mapping
        mapping_snd = [
                [x.shifted(shift[n]) for x in n_data]
                for n, n_data in enumerate(snd.id().mapping)
                ]
        inclusion_fst = OgMap(fst, disjoint_union, mapping_fst,
                              wfcheck=False)
        inclusion_snd = OgMap(snd, disjoint_union, mapping_snd,
                              wfcheck=False)

        return OgMapPair(inclusion_fst, inclusion_snd)

    @staticmethod
    def disjoint_union(fst, snd):
        """
        Returns the disjoint union of two OgPosets (the target
        of the coproduct cospan).
        """
        return OgPoset.coproduct(fst, snd).target

    @staticmethod
    def suspend(ogp, n=1):
        """ Returns the OgPoset suspended n times. """
        utils.typecheck(ogp, {'type': OgPoset})
        utils.typecheck(n, {
            'type': int,
            'st': lambda n: n >= 0,
            'why': 'expecting non-negative integer'})
        if n == 0:
            return ogp
        if n > 1:
            return OgPoset.suspend(
                    OgPoset.suspend(ogp, 1), n-1)
        face_data = [
                [
                    {'-': set(), '+': set()},
                    {'-': set(), '+': set()}
                ]] + ogp.face_data
        for x in ogp[0]:
            face_data[1][x.pos]['-'].add(0)
            face_data[1][x.pos]['+'].add(1)
        coface_data = [
                [
                    {'-': {x.pos for x in ogp[0]}, '+': set()},
                    {'-': set(), '+': {x.pos for x in ogp[0]}}
                ]] + ogp.coface_data
        return OgPoset(face_data, coface_data,
                       wfcheck=False, matchcheck=False)

    @staticmethod
    def gray(*ogps):
        """
        Returns the Gray product of a number of oriented graded posets.
        """
        if len(ogps) == 0:
            return OgPoset.point()
        if len(ogps) == 1:
            utils.typecheck(ogps[0], {'type': OgPoset})
            return ogps[0]
        if len(ogps) > 2:
            others = ogps[2:]
            return OgPoset.gray(OgPoset.gray(ogps[0], ogps[1]), *others)

        fst, snd = ogps[0], ogps[1]
        for x in (fst, snd):
            utils.typecheck(x, {'type': OgPoset})
        if len(fst) == 0 or len(snd) == 1:
            return fst
        if len(fst) == 1 or len(snd) == 0:
            return snd

        size1 = fst.size + [0 for _ in range(snd.dim)]
        size2 = snd.size + [0 for _ in range(fst.dim)]

        def pair(x, y):
            dim = x.dim + y.dim
            pos = y.pos + x.pos*size2[y.dim] + sum(
                    [size1[k]*size2[dim-k] for k in range(x.dim)])
            return El(dim, pos)

        def sndsign(n, sign):
            if n % 2 == 1:
                return utils.flip(sign)
            return sign

        face_data = [[] for _ in range(fst.dim + snd.dim + 1)]
        coface_data = [[] for _ in range(fst.dim + snd.dim + 1)]
        for x in fst:
            for y in snd:
                dim = x.dim + y.dim
                face_data[dim].append(
                    {sign:
                        {pair(z, y).pos
                         for z in fst.faces(x, sign)
                         }.union(
                            {pair(x, z).pos
                             for z in snd.faces(y, sndsign(x.dim, sign))
                             })
                     for sign in ('-', '+')})
                coface_data[dim].append(
                    {sign:
                        {pair(z, y).pos
                         for z in fst.cofaces(x, sign)
                         }.union(
                             {pair(x, z).pos
                              for z in snd.cofaces(y, sndsign(x.dim, sign))
                              })
                     for sign in ('-', '+')})
        return OgPoset(face_data, coface_data,
                       wfcheck=False, matchcheck=False)

    @staticmethod
    def join(*ogps):
        """ Returns the join of a number of oriented graded posets. """
        if len(ogps) == 0:
            return OgPoset.empty()
        if len(ogps) == 1:
            utils.typecheck(ogps[0], {'type': OgPoset})
            return ogps[0]
        if len(ogps) > 2:
            others = ogps[2:]
            return OgPoset.join(OgPoset.join(ogps[0], ogps[1]), *others)

        fst, snd = ogps[0], ogps[1]
        for x in (fst, snd):
            utils.typecheck(x, {'type': OgPoset})

        if len(fst) == 0:
            return snd
        if len(snd) == 0:
            return fst

        fst_bot_face_data = [[{'-': set(), '+': set()}]] + fst.face_data
        snd_bot_face_data = [[{'-': set(), '+': set()}]] + snd.face_data
        for x in fst_bot_face_data[1]:
            x['+'].add(0)
        for x in snd_bot_face_data[1]:
            x['+'].add(0)
        fst_bot_coface_data = [[
            {'-': set(), '+': {k for k in range(fst.size[0])}}
            ]] + fst.coface_data
        snd_bot_coface_data = [[
            {'-': set(), '+': {k for k in range(snd.size[0])}}
            ]] + snd.coface_data
        fst_bot = OgPoset(fst_bot_face_data, fst_bot_coface_data,
                          wfcheck=False, matchcheck=False)
        snd_bot = OgPoset(snd_bot_face_data, snd_bot_coface_data,
                          wfcheck=False, matchcheck=False)

        join_bot = OgPoset.gray(fst_bot, snd_bot)
        face_data = join_bot.face_data[1:]
        for x in face_data[0]:
            x['+'].clear()
        coface_data = join_bot.coface_data[1:]

        return OgPoset(face_data, coface_data,
                       wfcheck=False, matchcheck=False)

    # Private methods
    @staticmethod
    def __wfcheck(face_data):
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
            'why': 'expecting non-negative integer'})

        sizes = [len(_) for _ in face_data]
        for n, n_data in enumerate(face_data):
            # Check that faces are within bounds.
            k = max([i for x in n_data for sign in x for i in x[sign]],
                    default=-1)
            if (n == 0 and k >= 0) or k >= sizes[n-1]:
                raise ValueError(utils.value_err(k, 'out of bounds'))

            # Check that input/output are inhabited and disjoint.
            for i, x in enumerate(n_data):
                if not x['-'].isdisjoint(x['+']):
                    raise ValueError(utils.value_err(
                        face_data,
                        'input and output faces of El({}, {}) '
                        'are not disjoint'.format(
                            repr(n), repr(i))))
                if n > 0 and x['-'] == x['+'] == set():
                    raise ValueError(utils.value_err(
                        face_data,
                        'El({}, {}) must have at least one face'.format(
                            repr(n), repr(i))))

    @staticmethod
    def __coface_from_face(face_data):
        """
        Internal method constructing coface data from face data.
        Face data is presumed to be well-formed.
        """
        coface_data = [
                [
                    {'-': set(), '+': set()} for _ in n_data
                ]
                for n_data in face_data]
        for n, sn_data in enumerate(face_data[1:]):
            for k, x in enumerate(sn_data):
                for sign in '-', '+':
                    for i in x[sign]:
                        coface_data[n][i][sign].add(k)
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
        return "{}({})".format(
                type(self).__name__,
                ', '.join([repr(x) for x in self.as_list]))

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
        return type(self) is type(other) and \
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
            raise ValueError(utils.value_err(
                element, 'not in graded set'))

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

    def __init__(self, support, ambient,
                 wfcheck=True):
        if wfcheck:
            GrSubset.__wfcheck(support, ambient)

        self._support = support
        self._ambient = ambient

    def __eq__(self, other):
        return type(self) is type(other) and \
                self.support == other.support and \
                self.ambient == other.ambient

    def __str__(self):
        return '{} with {} elements in {}'.format(
                type(self).__name__,
                str(len(self.support)),
                str(self.ambient))

    def __contains__(self, item):
        return item in self.support

    def __len__(self):
        return len(self.support)

    def __getitem__(self, key):
        return GrSubset(self.support[key], self.ambient,
                        wfcheck=False)

    def __iter__(self):
        return iter(self.support)

    @property
    def support(self):
        """
        Returns the underlying GrSet (the 'support' of the subset).
        """
        return self._support

    @property
    def ambient(self):
        """ Returns the ambient OgPoset. """
        return self._ambient

    @property
    def dim(self):
        return self.support.dim

    @property
    def isclosed(self):
        """ Returns whether the subset is closed. """
        for n in range(self.dim, 0, -1):
            for x in self[n]:
                for face in self.ambient.faces(x):
                    if face not in self:
                        return False
        return True

    def union(self, *others):
        """
        Returns the union with other subsets of the same OgPoset.
        """
        others_support = []
        same_type = True
        for x in others:
            utils.typecheck(x, {
                'type': GrSubset,
                'st': lambda x: x.ambient == self.ambient,
                'why': 'not a subset of the same OgPoset'
                })
            if type(x) is not type(self):
                same_type = False
            others_support.append(x.support)

        union = self.support.union(*others_support)

        if same_type:  # return a Closed iff all are Closed
            return self.__class__(union, self.ambient,
                                  wfcheck=False)
        return GrSubset(union, self.ambient,
                        wfcheck=False)

    def intersection(self, *others):
        """
        Returns the intersection with other subsets of the same OgPoset.
        """
        others_support = []
        same_type = True
        for x in others:
            utils.typecheck(x, {
                'type': GrSubset,
                'st': lambda x: x.ambient == self.ambient,
                'why': 'not a subset of the same OgPoset'
                })
            if type(x) is not type(self):
                same_type = False
            others_support.append(x.support)

        intersection = self.support.intersection(*others_support)

        if same_type:  # return a Closed iff all are Closed
            return self.__class__(intersection, self.ambient,
                                  wfcheck=False)
        return GrSubset(intersection, self.ambient,
                        wfcheck=False)

    def issubset(self, other):
        utils.typecheck(other, {
            'type': GrSubset,
            'st': lambda x: x.ambient == self.ambient,
            'why': 'not a subset of the same OgPoset'
            })
        return self.support.issubset(other.support)

    def closure(self):
        """
        Returns the closure of the subset as an object of type Closed.
        """
        closure = self.support.copy()

        for n in range(self.dim, 0, -1):
            for element in closure[n]:
                for face in self.ambient.faces(element):
                    closure.add(face)

        return Closed(closure, self.ambient,
                      wfcheck=False)

    def image(self, ogmap):
        """
        Returns the image of the graded subset through an OgMap.
        """
        utils.typecheck(ogmap, {
            'type': OgMap,
            'st': lambda x: x.source == self.ambient,
            'why': 'OgMap source does not match ambient OgPoset'})

        image = GrSet()
        for x in self:
            if ogmap.isdefined(x):
                image.add(ogmap[x])

        # image of a Closed through an OgMap is Closed
        return self.__class__(image, ogmap.target,
                              wfcheck=False)

    # Internal methods
    @staticmethod
    def __wfcheck(support, ambient):
        utils.typecheck(support, {'type': GrSet})
        utils.typecheck(ambient, {'type': OgPoset})

        if support.dim > ambient.dim:
            raise ValueError(utils.value_err(
                support, 'does not define a subset'))
        for n in support.grades:
            if max([x.pos for x in support[n]]) >= ambient.size[n]:
                raise ValueError(utils.value_err(
                    support, 'does not define a subset'))


class Closed(GrSubset):
    """
    Subclass of GrSubset for closed subsets of oriented graded posets.
    """
    def __init__(self, support, ambient,
                 wfcheck=True):
        super().__init__(support, ambient, wfcheck)

        if wfcheck:
            if not self.isclosed:
                raise ValueError(utils.value_err(
                    support, 'not a closed subset'))

    @property
    def as_map(self):
        """
        Returns an injective OgMap representing the inclusion of
        the closed subset in the ambient.
        """
        mapping = [self.support[n].as_list
                   for n in range(self.dim + 1)]

        face_data = [
                [
                    {sign:
                        {mapping[n-1].index(y)
                         for y in self.ambient.faces(x, sign)}
                     for sign in ('-', '+')}
                    for x in n_data]
                for n, n_data in enumerate(mapping)]
        source = OgPoset.from_face_data(face_data, wfcheck=False)

        return OgMap(source, self.ambient, mapping,
                     wfcheck=False)

    def maximal(self):
        """
        Returns the subset of elements that are not below any other
        element in the graded set.
        """
        maximal = GrSet()
        for x in self:
            if self.ambient.cofaces(x).isdisjoint(
                    self.support[x.dim + 1]):
                maximal.add(x)
        return GrSubset(maximal, self.ambient,
                        wfcheck=False)

    def boundary_max(self, sign=None, dim=None):
        """
        Returns the set of maximal elements of the n-boundary of
        the closed set.
        """
        _sign = utils.flip(
                utils.mksign(sign)) if sign is not None else '-'
        dim = self.dim - 1 if dim is None else dim

        boundary_max = self.maximal().support[:dim]

        for x in self[dim]:
            if self.ambient.cofaces(x, _sign).isdisjoint(
                    self.support[x.dim + 1]):
                boundary_max.add(x)
            if sign is None and self.ambient.cofaces(x, '+').isdisjoint(
                    self.support[x.dim + 1]):
                boundary_max.add(x)

        return GrSubset(boundary_max, self.ambient,
                        wfcheck=False)

    def boundary(self, sign=None, dim=None):
        """
        Returns the n-boundary of the closed set.
        """
        if isinstance(dim, int) and dim >= self.dim:
            return self
        return self.boundary_max(sign, dim).closure()

    @staticmethod
    def subset(grsubset,
               wfcheck=True):
        if wfcheck:
            if not grsubset.isclosed:
                raise ValueError(grsubset.support, 'not a closed subset')
        return Closed(grsubset.support, grsubset.ambient,
                      wfcheck=False)


class OgMap:
    """
    Class for (partial) maps of oriented graded posets.
    """

    def __init__(self, source, target, mapping=None,
                 wfcheck=True):
        if wfcheck:
            OgMap.__wfcheck(source, target, mapping)

        self._source = source
        self._target = target
        if mapping is None:
            mapping = [[None for _ in range(source.size[n])]
                       for n in range(len(source.size))]
        self._mapping = mapping

    def __eq__(self, other):
        return type(self) == type(other) and \
                self.source == other.source and \
                self.target == other.target and \
                self.mapping == other.mapping

    def __str__(self):
        return '{} from {} to {}'.format(
                type(self).__name__,
                str(self.source), str(self.target))

    def __getitem__(self, element):
        if element in self.source:
            return self.mapping[element.dim][element.pos]
        raise ValueError(utils.value_err(
            element, 'not in source'))

    def __setitem__(self, element, image):
        if element in self.source:
            self.__extensioncheck(element, image)
            self._mapping[element.dim][element.pos] = image
        else:
            raise ValueError(utils.value_err(
                element, 'not in source'))

    @property
    def source(self):
        return self._source

    @property
    def target(self):
        return self._target

    @property
    def mapping(self):
        return self._mapping

    @property
    def istotal(self):
        """ Returns whether the map is total. """
        for n_data in self.mapping:
            if not all(n_data):
                return False
        return True

    @property
    def isinjective(self):
        """ Returns whether the map is injective. """
        image_list = [x for n_data in self.mapping for x in n_data
                      if x is not None]
        if len(image_list) == len(set(image_list)):
            return True
        return False

    @property
    def issurjective(self):
        """ Returns whether the map is surjective. """
        image_set = {x for n_data in self.mapping for x in n_data
                     if x is not None}
        if len(image_set) == len(self.target):
            return True
        return False

    @property
    def isiso(self):
        """ Returns whether the map is an isomorphism """
        return self.istotal and self.isinjective and self.issurjective

    def isdefined(self, element):
        """ Returns whether the map is defined on an element. """
        if element in self.source and self[element] is not None:
            return True
        return False

    def then(self, other, *others):
        """ Returns the composite with other maps. """
        if len(others) > 0:
            return self.then(other).then(*others)

        if isinstance(other, OgMapPair):
            return OgMapPair(
                    self.then(other.fst),
                    self.then(other.snd))

        utils.typecheck(other, {
            'type': OgMap,
            'st': lambda x: x.source == self.target,
            'why': 'source does not match target of first map'})
        mapping = [
                [other.mapping[x.dim][x.pos] if x is not None
                 else None for x in n_data]
                for n_data in self.mapping]

        if type(self) is type(other):  # for inheritance
            return self.__class__(self.source, other.target, mapping,
                                  wfcheck=False)
        return OgMap(self.source, other.target, mapping,
                     wfcheck=False)

    def inv(self):
        """ Returns the inverse of the OgMap if it is an isomorphism. """
        if not self.isiso:
            raise ValueError(utils.value_err(
                self, 'not an isomorphism'))
        mapping_inv = [[None for _ in n_data]
                       for n_data in self.mapping]

        for x in self.source:
            image = self[x]
            mapping_inv[image.dim][image.pos] = x

        return OgMap(self.target, self.source, mapping_inv,
                     wfcheck=False)

    @staticmethod
    def compose(one, *others):
        utils.typecheck(one, {'type': OgMap})
        if len(others) > 0:
            return one.then(*others)
        return one

    # Private methods.
    def __extensioncheck(self, element, image):
        if image not in self.target:
            raise ValueError(utils.value_err(
                image, 'not in target'))
        if self.isdefined(element):
            raise ValueError(utils.value_err(
                element, 'already defined on element'))
        if image.dim > element.dim:
            raise ValueError(utils.value_err(
                image, 'exceeds dimension of {}'.format(
                    repr(element))))

        el_underset = GrSubset(GrSet(element), self.source).closure()

        for x in el_underset[:element.dim]:
            if not self.isdefined(x):
                raise ValueError(utils.value_err(
                    element, 'map undefined on {} below {}'.format(
                        repr(x), repr(element))))

        img_underset = GrSubset(GrSet(image), self.target).closure()
        for n in range(element.dim):
            for sign in '-', '+':
                if el_underset.boundary(sign, n).image(self) != \
                        img_underset.boundary(sign, n):
                    raise ValueError(utils.value_err(
                        image,
                        'assignment does not respect ({}, {})-boundary'.format(
                            sign, repr(n))))

    @staticmethod
    def __wfcheck(source, target, mapping):
        for x in source, target:
            utils.typecheck(x, {'type': OgPoset})
        if mapping is not None:  # otherwise nothing else to check
            utils.typecheck(mapping, {'type': list}, {'type': list})

            mapping_size = [len(_) for _ in mapping]
            if mapping_size != source.size:
                raise ValueError(utils.value_err(
                    mapping, 'wrong size'))

            check_map = OgMap(source, target)
            # Extend check_map one element at a time according to data in
            # mapping, if this gives no error the check is passed.
            for x in source:
                if mapping[x.dim][x.pos] is not None:
                    check_map[x] = mapping[x.dim][x.pos]


class OgMapPair(tuple):
    """
    Class for pairs of maps.
    """

    def __new__(self, fst, snd):
        for x in fst, snd:
            utils.typecheck(x, {'type': OgMap})
        return tuple.__new__(OgMapPair, (fst, snd))

    def __str__(self):
        return '({}, {})'.format(str(self.fst), str(self.snd))

    def __eq__(self, other):
        return type(self) is type(other) and \
                self.fst == other.fst and self.snd == other.snd

    @property
    def fst(self):
        return self[0]

    @property
    def snd(self):
        return self[1]

    @property
    def source(self):
        if self.isspan:
            return self.fst.source
        return self.fst.source, self.snd.source

    @property
    def target(self):
        if self.iscospan:
            return self.fst.target
        return self.fst.target, self.snd.target

    @property
    def isspan(self):
        return self.fst.source == self.snd.source

    @property
    def iscospan(self):
        return self.fst.target == self.snd.target

    @property
    def isparallel(self):
        return self.isspan and self.iscospan

    @property
    def istotal(self):
        return self.fst.istotal and self.snd.istotal

    @property
    def isinjective(self):
        return self.fst.isinjective and self.snd.isinjective

    def then(self, other, *others):
        """ Returns the composite with other pairs or maps. """
        if len(others) > 0:
            return self.then(other).then(*others)

        if isinstance(other, OgMapPair):
            return OgMapPair(
                    self.fst.then(other.fst),
                    self.snd.then(other.snd))

        return OgMapPair(
                self.fst.then(other),
                self.snd.then(other))

    def coequaliser(self,
                    wfcheck=True):
        """
        Returns the coequaliser of a pair of injective total maps,
        if it exists.
        """
        if not (self.isparallel and self.istotal and self.isinjective):
            raise ValueError(utils.value_err(
                self,
                'expecting a parallel pair of injective total maps'))

        mapping = self.target.id().mapping
        to_delete = GrSet()
        for x in self.source:
            if self.fst[x] != self.snd[x]:
                to_delete.add(self.snd[x])
                mapping[self.snd[x].dim][self.snd[x].pos] = self.fst[x]

        # Shift assignments following deletion of elements
        shift_list = [[0 for _ in n_data]
                      for n_data in mapping]
        for x in to_delete:
            for k, y in enumerate(mapping[x.dim]):
                if y.pos >= x.pos:
                    shift_list[x.dim][k] -= 1
        mapping = [
                [x.shifted(k) for x, k in zip(n_data, n_shifts)]
                for n_data, n_shifts in zip(mapping, shift_list)]

        face_data = [
                [
                    {sign:
                        {mapping[n-1][j].pos for j in x[sign]}
                     for sign in x}
                    for k, x in enumerate(n_data) if El(n, k) not in to_delete
                ]
                for n, n_data in enumerate(self.target.face_data)]
        quotient = OgPoset.from_face_data(face_data, wfcheck)

        return OgMap(self.target, quotient, mapping, wfcheck=False)

    def pushout(self,
                wfcheck=True):
        """
        Returns the pushout of a span of injective total maps, if it
        exists.
        """
        if not (self.isspan and self.istotal and self.isinjective):
            raise ValueError(utils.value_err(
                self,
                'expecting a span of injective total maps'))

        coproduct = OgPoset.coproduct(self.fst.target, self.snd.target)
        coequaliser = self.then(coproduct).coequaliser(wfcheck)

        return coproduct.then(coequaliser)
