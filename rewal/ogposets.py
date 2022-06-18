"""
Implements oriented graded posets, their elements, subsets, and maps.
"""

import numpy as np

from rewal import utils


class El(tuple):
    """
    Class for elements of an oriented graded poset, encoded as pairs
    of non-negative integers: the dimension of the element, and its
    position in a linear order of elements of the given dimension.

    Arguments
    ---------
    dim : int
        The dimension of the element.
    pos : int
        The position of the element.

    Examples
    --------
    >>> x = El(2, 3)
    >>> assert x.dim == 2 and x.pos == 3
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
        return isinstance(other, El) and \
                self.dim == other.dim and self.pos == other.pos

    def __hash__(self):
        return hash(repr(self))

    @property
    def dim(self):
        """
        Returns the dimension of the element.

        Returns
        -------
        dim : int
            The dimension of the element.
        """
        return self[0]

    @property
    def pos(self):
        """
        Returns the position of the element.

        Returns
        -------
        pos : int
            The position of the element
        """
        return self[1]

    def shifted(self, k):
        """
        Returns the element of the same dimension, with position shifted
        by a given integer.

        Parameters
        ----------
        k : int
            The shift in position.

        Returns
        -------
        el : :class:`El`
            The shifted element.
        """
        utils.typecheck(k, {
            'type': int,
            'st': lambda n: self.pos + n >= 0,
            'why': 'shifted position must be non-negative'
            })
        return El(self.dim, self.pos + k)


class OgPoset:
    """
    Class for oriented graded posets, that is, finite graded posets
    with an orientation, defined as a :code:`{'-', '+'}`-labelling of the
    edges of their Hasse diagram.

    In this implementation, the elements of a given dimension (grade)
    are linearly ordered, so that each element is identified by its
    dimension and the position in the linear order, encoded as an object
    of class :class:`El`.

    If :code:`El(n, k)` covers :code:`El(n-1, j)` with orientation
    :code:`o`, we say that :code:`El(n-1, j)` is an *input face* of
    :code:`El(n, k)` if :code:`o == '-'` and an *output face* of
    :code:`El(n, k)` if :code:`o == '+'`.

    Defining an :class:`OgPoset` directly is not recommended; use
    constructors of :class:`shapes.Shape` instead.

    Arguments
    ---------
    face_data : list of list of dict of set of int
        Data encoding the oriented graded poset as follows:
        :code:`j in face_data[n][k][o]` if and only if
        :code:`El(n, k)` covers :code:`El(n-1, j)` with orientation
        :code:`o`, where :code:`o == '-'` or :code:`o == '+'`.

    coface_data: list of list of dict of set of int
        Data encoding the oriented graded poset as follows:
        :code:`j in coface_data[n][k][o]` if and only if
        :code:`El(n+1, j)` covers :code:`El(n, k)` with orientation
        :code:`o`, where :code:`o == '-'` or :code:`o == '+'`.

    Keyword arguments
    -----------------
    wfcheck : bool
        Check that the data is well-formed (default is :code:`True`)
    matchcheck : bool
        Check that `face_data` and `coface_data` match
        (default is :code:`True`)

    Notes
    -----
    Each of `face_data`, `coface_data` determines the other uniquely.
    There is an alternative constructor :meth:`from_face_data` that computes
    `coface_data` from `face_data`.

    Examples
    --------
    Let us construct explicitly the “oriented face poset” of an arrow, or
    directed edge.

    >>> face_data = [
    ...    [
    ...        {'-': set(), '+': set()},
    ...        {'-': set(), '+': set()},
    ...    ], [
    ...        {'-': {0}, '+': {1}}
    ...    ]]
    >>> coface_data = [
    ...    [
    ...        {'-': {0}, '+': set()},
    ...        {'-': set(), '+': {0}},
    ...    ], [
    ...        {'-': set(), '+': set()}
    ...    ]]
    >>> arrow = OgPoset(face_data, coface_data)

    This has two 0-dimensional elements and one 1-dimensional element.

    >>> arrow.size
    [2, 1]

    We can visualise its Hasse diagram, with orientation conveyed by colour
    (magenta for input, blue for output) and direction of arrows.

    >>> arrow.hasse(path='docs/_static/img/OgPoset_arrow.png')

    .. image:: ../_static/img/OgPoset_arrow.png
        :align: center

    We can ask for the faces and cofaces of a specific element.

    >>> arrow.faces(El(1, 0), '-')
    GrSet(El(0, 0))
    >>> arrow.cofaces(El(0, 1))
    GrSet(El(1, 0))

    We can construct other oriented graded posets using various operations,
    such as suspensions, Gray products, joins, or duals.

    >>> print(arrow.suspend())
    OgPoset with [2, 2, 1] elements
    >>> print(arrow * arrow)
    OgPoset with [4, 4, 1] elements
    >>> print(arrow >> arrow)
    OgPoset with [4, 6, 4, 1] elements
    >>> print(arrow.dual())
    OgPoset with [2, 1] elements
    """

    def __init__(self, face_data, coface_data, **params):
        wfcheck = params.get('wfcheck', True)
        matchcheck = params.get('matchcheck', True)
        if wfcheck:
            OgPoset._wfcheck(face_data)

        if matchcheck:
            if not coface_data == OgPoset._coface_from_face(face_data):
                raise ValueError(utils.value_err(
                    coface_data, 'face and coface data do not match'))

        self._face_data = face_data
        self._coface_data = coface_data

        # Enable method chaining syntax
        self.gray = self._gray
        self.join = self._join
        self.suspend = self._suspend
        self.dual = self._dual
        if hasattr(self, 'atom'):
            self.atom = self._atom
        if hasattr(self, 'paste'):
            self.paste = self._paste

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
        return isinstance(other, OgPoset) and \
                self.face_data == other.face_data and \
                self.coface_data == other.coface_data

    def __add__(self, other):
        return OgPoset.disjoint_union(self, other)

    def __mul__(self, other):
        """ Returns the Gray product. """
        return self.gray(other)

    def __pow__(self, times):
        utils.typecheck(times, {'type': int})
        return self.__class__.gray(*[self for _ in range(times)])

    def __rshift__(self, other):
        """ Returns the join. """
        return self.join(other)

    def __lshift__(self, other):
        return other.join(self)

    def __invert__(self):
        """ Returns the dual. """
        return self.dual()

    @property
    def face_data(self):
        """
        Returns the face data as given to the object constructor.

        An :class:`OgPoset` is meant to be immutable; create a new
        object if you need to modify the face data.

        Returns
        -------
        face_data : list of list of dict of set of int
            The face data as given to the object constructor.
        """
        return self._face_data

    @property
    def coface_data(self):
        """
        Returns the coface data as given to the object constructor.

        An :class:`OgPoset` is meant to be immutable; create a new
        object if you need to modify the coface data.

        Returns
        -------
        coface_data : list of list of dict of set of int
            The coface data as given to the object constructor.
        """
        return self._coface_data

    @property
    def size(self):
        """
        Returns the number of elements in each dimension as a list.

        Returns
        -------
        size : list of int
            The :code:`k` th entry is the number of
            :code:`k` -dimensional elements.
        """
        return [len(_) for _ in self.face_data]

    @property
    def dim(self):
        """
        Returns the dimension of the object, that is, the maximum of
        the dimensions of its elements.

        Returns
        -------
        dim : int
            The dimension of the object.
        """
        return len(self.face_data) - 1

    @property
    def as_chain(self):
        """
        Returns a “chain complex” representation of the face data.

        Returns
        -------
        chain : list of :class:`numpy.array`
            Encodes the face data as follows:
            :code:`chain[n][i][j] == 1` if
            :code:`El(n, i)` is an output face of :code:`El(n+1, j)`,
            :code:`-1` if it is an input face, :code:`0` otherwise.
        """
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
        """
        Returns the closed subset of all elements.

        Returns
        -------
        all : :class:`Closed`
            The closed subset of all elements of the object.
        """
        return Closed(
                GrSet(*[El(n, k) for n, n_data in enumerate(self.face_data)
                        for k, x in enumerate(n_data)]),
                self,
                wfcheck=False)

    def none(self):
        """
        Returns the empty closed subset.

        Returns
        -------
        none : :class:`Closed`
            The closed subset with no elements.
        """
        return Closed(GrSet(), self,
                      wfcheck=False)

    def underset(self, *elements):
        """
        Returns the closure of a set of elements in the object.

        Arguments
        ---------
        elements : :class:`El`
            Any number of elements.

        Returns
        -------
        underset : :class:`Closed`
            The downwards closure of the given elements.
        """
        return GrSubset(GrSet(*elements), self).closure()

    def maximal(self):
        """
        Returns the subset of maximal elements, that is, those
        that are not covered by any elements.

        Returns
        -------
        maximal : :class:`GrSubset`
            The subset of maximal elements.
        """
        maximal = GrSet()
        for x in self:
            if self.cofaces(x) == GrSet():
                maximal.add(x)
        return GrSubset(maximal, self,
                        wfcheck=False)

    def faces(self, element, sign=None):
        """
        Returns the faces of an element as a graded set.

        Arguments
        ---------
        element : :class:`El`
            An element of the object.
        sign : str or int, optional
            Orientation: :code:`'-'` for input, :code:`'+'` for output,
            :code:`None` (default) for both.

        Returns
        -------
        faces : :class:`GrSet`
            The set of faces of the given element.
        """
        if sign is None:
            return self.faces(element, '-').union(
                self.faces(element, '+'))
        sign = utils.mksign(sign)
        utils.typecheck(element, {
            'type': El,
            'st': lambda x: x.dim <= self.dim and x.pos <= self.size[x.dim],
            'why': 'out of bounds'})
        if element.dim == 0:
            return GrSet()
        return GrSet(
                *[El(element.dim-1, i)
                  for i in self.face_data[element.dim][element.pos][sign]]
                )

    def cofaces(self, element, sign=None):
        """
        Returns the cofaces of an element as a graded set.

        Arguments
        ---------
        element : :class:`El`
            An element of the object.
        sign : str or int, optional
            Orientation: :code:`'-'` for input, :code:`'+'` for output,
            :code:`None` (default) for both.

        Returns
        -------
        cofaces : :class:`GrSet`
            The set of cofaces of the given element.
        """
        if sign is None:
            return self.cofaces(element, '-').union(
                self.cofaces(element, '+'))
        sign = utils.mksign(sign)
        utils.typecheck(element, {
            'type': El,
            'st': lambda x: x.dim <= self.dim and x.pos <= self.size[x.dim],
            'why': 'out of bounds'})
        if element.dim == self.dim:
            return GrSet()
        return GrSet(
                *[El(element.dim + 1, i)
                  for i in self.coface_data[element.dim][element.pos][sign]]
                )

    def id(self):
        """
        Returns the identity map on the object.

        Returns
        -------
        id : :class:`OgMap`
            The identity map on the object.
        """
        mapping = [
                [
                    El(n, k) for k, x in enumerate(n_data)
                ]
                for n, n_data in enumerate(self.face_data)
                ]
        return OgMap(self, self, mapping,
                     wfcheck=False)

    def image(self, ogmap):
        """
        Returns the image of the object through a map.

        Arguments
        ---------
        ogmap : :class:`OgMap`
            A map from the object to another :class:`OgPoset`.

        Returns
        -------
        image : :class:`Closed`
            The image of the object through the given map.
        """
        return self.all().image(ogmap)

    def boundary(self, sign=None, dim=None):
        """
        Returns the inclusion of the boundary of a given orientation
        and dimension into the object.

        Arguments
        ---------
        sign : str or int
            Orientation: :code:`'-'` for input, :code:`'+'` for output,
            :code:`None` (default) for both.
        dim : int
            Dimension of the boundary (default is :code:`self.dim - 1`).

        Returns
        -------
        boundary : :class:`OgMap`
            The inclusion of the requested boundary into the object.
        """
        if isinstance(dim, int) and dim >= self.dim:
            return self.id()
        return self.all().boundary(sign, dim).as_map

    @property
    def input(self):
        """
        Alias for :code:`boundary('-')`.
        """
        return self.boundary('-')

    @property
    def output(self):
        """
        Alias for :code:`boundary('+')`.
        """
        return self.boundary('+')

    @classmethod
    def from_face_data(cls, face_data,
                       wfcheck=True):
        """
        Alternative constructor computing `coface_data` from `face_data`.

        Arguments
        ---------
        face_data : list of list of dict of set of int
            As in the main constructor.
        wfcheck : bool
            Check that the data is well-formed (default is :code:`True`).
        """
        if wfcheck:
            cls._wfcheck(face_data)
        coface_data = cls._coface_from_face(face_data)

        return cls(face_data, coface_data,
                   wfcheck=False, matchcheck=False)

    @staticmethod
    def empty():
        """
        Returns the initial oriented graded poset, with no elements.

        Returns
        -------
        empty : :class:`OgPoset`
            The empty oriented graded poset.
        """
        return OgPoset([], [],
                       wfcheck=False, matchcheck=False)

    @staticmethod
    def point():
        """
        Returns the terminal oriented graded poset, with a single element.

        Returns
        -------
        point : :class:`OgPoset`
            The oriented graded poset with a single element.
        """
        return OgPoset(
            [[{'-': set(), '+': set()}]],
            [[{'-': set(), '+': set()}]],
            wfcheck=False, matchcheck=False)

    @staticmethod
    def coproduct(fst, snd):
        """
        Returns the coproduct cospan of two oriented graded posets.

        Arguments
        ---------
        fst : :class:`OgPoset`
            The first factor of the coproduct.
        snd : :class:`OgPoset`
            The second factor of the coproduct.

        Returns
        -------
        coproduct : :class:`OgMapPair`
            The coproduct cospan.
        """
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
        Returns the disjoint union of two oriented graded posets, that is,
        the target of their coproduct cospan.

        This method can be called with the math operator :code:`+`, that is,
        :code:`fst + snd` is equivalent to :code:`disjoint_union(fst, snd)`.

        Arguments
        ---------
        fst : :class:`OgPoset`
            The first factor of the disjoint union.
        snd : :class:`OgPoset`
            The second factor of the disjoint union.

        Returns
        -------
        disjoint_union : :class:`OgPoset`
            The disjoint union of the two.
        """
        return OgPoset.coproduct(fst, snd).target

    @staticmethod
    def suspend(ogp, n=1):
        """
        Returns the n-fold suspension of an oriented graded poset.

        This static method can be also used as a bound method, that is,
        :code:`self.suspend(n)` is equivalent to :code:`suspend(self, n)`.

        Arguments
        ---------
        ogp : :class:`OgPoset`
            The object to suspend.
        n : int
            The number of iterations of the suspension (default is 1).

        Returns
        -------
        suspension : :class:`OgPoset`
            The suspended object.
        """
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
                ], *ogp.face_data]
        for x in ogp[0]:
            face_data[1][x.pos]['-'].add(0)
            face_data[1][x.pos]['+'].add(1)
        coface_data = [
                [
                    {'-': {x.pos for x in ogp[0]}, '+': set()},
                    {'-': set(), '+': {x.pos for x in ogp[0]}}
                ], *ogp.coface_data]
        return OgPoset(face_data, coface_data,
                       wfcheck=False, matchcheck=False)

    def _suspend(self, n=1):
        return self.__class__.suspend(self, n)

    @staticmethod
    def gray(*ogps):
        """
        Returns the Gray product of any number of oriented graded posets.

        This method can be called with the math operator :code:`*`, that is,
        :code:`fst * snd` is equivalent to :code:`gray(fst, snd)`.

        This static method can also be used as a bound method, that is,
        :code:`fst.gray(*ogps)` is equivalent to :code:`gray(fst, *ogps)`.

        Arguments
        ---------
        ogps : :class:`OgPoset`
            Any number of oriented graded posets.

        Returns
        -------
        gray : :class:`OgPoset`
            The Gray product of the arguments.
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

    def _gray(self, *others):
        return self.__class__.gray(self, *others)

    def bot(self):
        """
        Returns the object augmented with a bottom element, covered
        with orientation :code:`'+'`.

        Returns
        -------
        bot : :class:`OgPoset`
            The object augmented with a bottom element.
        """
        if len(self) == 0:
            return OgPoset.point()
        face_data = [[{'-': set(), '+': set()}], *self.face_data]
        for x in face_data[1]:
            x['+'].add(0)
        coface_data = [[
            {'-': set(), '+': {k for k in range(self.size[0])}}
            ], *self.coface_data]
        return OgPoset(face_data, coface_data,
                       wfcheck=False, matchcheck=False)

    @staticmethod
    def join(*ogps):
        """
        Returns the join of any number of oriented graded posets.

        This method can be called with the bitwise operators :code:`>>`
        and :code:`<<`, that is, :code:`fst >> snd` is equivalent to
        :code:`join(fst, snd)` and :code:`fst << snd` is equivalent to
        :code:`join(snd, fst)`.

        This static method can also be used as a bound method, that is,
        :code:`fst.join(*ogps)` is equivalent to :code:`join(fst, *ogps)`.

        Arguments
        ---------
        ogps : :class:`OgPoset`
            Any number of oriented graded posets.

        Returns
        -------
        gray : :class:`OgPoset`
            The join of the arguments.
        """
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

        join_bot = OgPoset.gray(
                fst.bot(),
                snd.bot())
        face_data = join_bot.face_data[1:]
        for x in face_data[0]:
            x['+'].clear()
        coface_data = join_bot.coface_data[1:]

        return OgPoset(face_data, coface_data,
                       wfcheck=False, matchcheck=False)

    def _join(self, *others):
        return self.__class__.join(self, *others)

    @staticmethod
    def dual(ogp, *dims):
        """
        Returns an oriented graded poset with orientations reversed
        in given dimensions.

        This static method can be also used as a bound method, that is,
        :code:`self.dual(*dims)` is equivalent to :code:`dual(self, *dims)`.

        Arguments
        ---------
        ogp : :class:`OgPoset`
            An oriented graded poset.
        dims : int
            Any number of dimensions; if none, defaults to *all* dimensions.

        Returns
        -------
        dual : :class:`OgPoset`
            The oriented graded poset, dualised in the given dimensions.
        """
        utils.typecheck(ogp, {'type': OgPoset})
        for x in dims:
            utils.typecheck(x, {'type': int})
        if len(dims) == 0:  # default is dualising in all dimensions
            dims = range(ogp.dim + 1)

        face_data = [
                [
                    {sign: x[utils.flip(sign)] if n in dims else x[sign]
                     for sign in ('-', '+')}
                    for x in n_data
                ]
                for n, n_data in enumerate(ogp.face_data)]
        coface_data = [
                [
                    {sign: x[utils.flip(sign)] if n+1 in dims else x[sign]
                     for sign in ('-', '+')}
                    for x in n_data
                ]
                for n, n_data in enumerate(ogp.coface_data)]
        return OgPoset(face_data, coface_data,
                       wfcheck=False, matchcheck=True)

    def _dual(self, *dims):
        return self.__class__.dual(self, *dims)

    def op(self):
        """
        Returns the :meth:`dual` in all *odd* dimensions.
        """
        odds = [n for n in range(self.dim + 1) if n % 2 == 1]
        return self.dual(*odds)

    def co(self):
        """
        Returns the :meth:`dual` in all *even* dimensions.
        """
        evens = [n for n in range(self.dim + 1) if n % 2 == 0]
        return self.dual(*evens)

    def hasse(self, **params):
        """
        Alias for :code:`hasse.draw(self, **params)`.
        """
        from rewal.hasse import draw
        return draw(self, **params)

    # Private methods
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
            'why': 'expecting non-negative integer'})

        sizes = [len(n_data) for n_data in face_data]
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
    def _coface_from_face(face_data):
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
    Class for sets of elements of an oriented graded poset, graded
    by their dimension.
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

    def difference(self, other):
        """
        Returns the difference of the graded set with another graded set.
        """
        utils.typecheck(other, {'type': GrSet})

        difference_as_set = self.as_set.difference(other.as_set)
        return GrSet(*difference_as_set)

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

    def __init__(self, support, ambient, **params):
        wfcheck = params.get('wfcheck', True)
        if wfcheck:
            GrSubset._wfcheck(support, ambient)

        self._support = support
        self._ambient = ambient

    def __eq__(self, other):
        return isinstance(other, GrSubset) and \
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

    def difference(self, other):
        """
        Returns the difference with another subset of the same OgPoset.
        """
        utils.typecheck(other, {
            'type': GrSubset,
            'st': lambda x: x.ambient == self.ambient,
            'why': 'not a subset of the same OgPoset'
            })
        difference = self.support.difference(other.support)
        return GrSubset(difference, self.ambient,
                        wfcheck=False)

    def issubset(self, other):
        utils.typecheck(other, {
            'type': GrSubset,
            'st': lambda x: x.ambient == self.ambient,
            'why': 'not a subset of the same OgPoset'
            })
        return self.support.issubset(other.support)

    def isdisjoint(self, other):
        utils.typecheck(other, {
            'type': GrSubset,
            'st': lambda x: x.ambient == self.ambient,
            'why': 'not a subset of the same OgPoset'
            })
        return self.support.isdisjoint(other.support)

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
    def _wfcheck(support, ambient):
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
    def __init__(self, support, ambient, **params):
        wfcheck = params.get('wfcheck', True)
        super().__init__(support, ambient, **params)

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

    @property
    def ispure(self):
        """
        Returns whether the maximal elements of the closed subset
        all have the same dimension.
        """
        for x in self.maximal():
            if x.dim != self.dim:
                return False
        return True

    @property
    def isround(self):
        """
        Returns whether the closed subset "has spherical boundary".
        """
        if not self.ispure:
            return False
        boundary_in = self.boundary('-')
        boundary_out = self.boundary('+')
        intersection = boundary_in.intersection(boundary_out)
        for k in range(self.dim-2, -1, -1):
            boundary_in = boundary_in.boundary('-')
            boundary_out = boundary_out.boundary('+')
            if not intersection.issubset(boundary_in.union(boundary_out)):
                return False
            intersection = boundary_in.intersection(boundary_out)
        return True

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
        dim = self.dim-1 if dim is None else dim

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

    @property
    def input(self):
        return self.boundary('-')

    @property
    def output(self):
        return self.boundary('+')

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

    def __init__(self, source, target, mapping=None, **params):
        wfcheck = params.get('wfcheck', True)
        if wfcheck:
            OgMap._wfcheck(source, target, mapping)

        self._source = source
        self._target = target
        if mapping is None:
            mapping = [[None for _ in n_data]
                       for n_data in self.source.face_data]
        self._mapping = mapping

        # Enable method chaining syntax
        self.gray = self._gray
        self.join = self._join
        self.dual = self._dual

    def __eq__(self, other):
        return isinstance(other, OgMap) and \
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
        if element not in self.source:
            raise ValueError(utils.value_err(
                element, 'not in source'))
        self._extensioncheck(element, image)
        self._mapping[element.dim][element.pos] = image

    def __call__(self, other):
        if isinstance(other, El):
            return self[other]
        if isinstance(other, GrSubset):
            return other.image(self)
        if isinstance(other, GrSet):
            subset = GrSubset(other, self.source)
            return subset.image(self)
        raise TypeError(utils.type_err(
            El, other))

    def __mul__(self, other):
        return self.gray(other)

    def __pow__(self, times):
        utils.typecheck(times, {'type': int})
        return self.__class__.gray(*[self for _ in range(times)])

    def __rshift__(self, other):
        return self.join(other)

    def __lshift__(self, other):
        return other.join(self)

    def __invert__(self):
        return self.dual()

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

    def image(self):
        """ The image of the source through the mapping. """
        return self.source.all().image(self)

    def boundary(self, sign=None, dim=None):
        """
        The map restricted to a boundary of its source.
        """
        return self.source.boundary(
                sign, dim).then(self)

    @property
    def input(self):
        return self.boundary('-')

    @property
    def output(self):
        return self.boundary('+')

    def bot(self):
        """
        Extension of OgPoset.bot to maps.
        """
        source = self.source.bot()
        target = self.target.bot()
        mapping = [[El(0, 0)]] + [
                [El(x.dim + 1, x.pos) for x in n_data]
                for n_data in self.mapping]

        return OgMap(source, target, mapping,
                     wfcheck=False)

    @staticmethod
    def gray(*maps):
        for f in maps:
            utils.typecheck(f, {'type': OgMap})
        if len(maps) == 0:
            return OgMap.point().id()
        if len(maps) == 1:
            return maps[0]

        fst, snd = maps[0], maps[1]
        if len(maps) > 2:
            return OgMap.gray(
                    OgMap.gray(fst, snd), *maps[2:])

        size1 = fst.target.size + [0 for _ in range(snd.target.dim)]
        size2 = snd.target.size + [0 for _ in range(fst.target.dim)]

        def pair(x, y):
            if x is None or y is None:
                return None
            dim = x.dim + y.dim
            pos = y.pos + x.pos*size2[y.dim] + sum(
                    [size1[k]*size2[dim-k] for k in range(x.dim)])
            return El(dim, pos)

        mapping = [[] for _ in range(fst.source.dim + snd.source.dim + 1)]
        for x in fst.source:
            for y in snd.source:
                mapping[x.dim + y.dim].append(
                    pair(fst[x], snd[y]))

        return OgMap(
                OgPoset.gray(fst.source, snd.source),
                OgPoset.gray(fst.target, snd.target),
                mapping,
                wfcheck=False)

    def _gray(self, *others):
        return self.__class__.gray(self, *others)

    @staticmethod
    def join(*maps):
        for f in maps:
            utils.typecheck(f, {'type': OgMap})
        if len(maps) == 0:
            return OgPoset.empty().id()
        if len(maps) == 1:
            return maps[0]

        fst, snd = maps[0], maps[1]
        if len(maps) > 2:
            return OgMap.join(
                    OgMap.join(fst, snd), *maps[2:])

        join_bot = OgMap.gray(
                fst.bot(), snd.bot())
        mapping = [
                [El(x.dim-1, x.pos) for x in n_data]
                for n_data in join_bot.mapping[1:]]

        return OgMap(
                OgPoset.join(fst.source, snd.source),
                OgPoset.join(fst.target, snd.target),
                mapping,
                wfcheck=False)

    def _join(self, *others):
        return self.__class__.join(self, *others)

    @staticmethod
    def dual(ogmap, *dims):
        utils.typecheck(ogmap, {'type': OgMap})
        return OgMap(OgPoset.dual(ogmap.source, *dims),
                     OgPoset.dual(ogmap.target, *dims),
                     ogmap.mapping,
                     wfcheck=False)

    def _dual(self, *dims):
        return self.__class__.dual(self, *dims)

    def op(self):
        odds = [n for n in range(self.dim + 1) if n % 2 == 1]
        return self.dual(*odds)

    def co(self):
        evens = [n for n in range(self.dim + 1) if n % 2 == 0]
        return self.dual(*evens)

    def hasse(self, **params):
        from rewal.hasse import draw
        return draw(self, **params)

    # Private methods.
    def _extensioncheck(self, element, image):
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

        el_underset = self.source.underset(element)

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
                        'assignment does not respect '
                        '({}, {})-boundary of {}'.format(
                            sign, repr(n), repr(element))))

    @staticmethod
    def _wfcheck(source, target, mapping):
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
        return isinstance(other, OgMapPair) and \
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
        Returns the coequaliser of a pair of total maps, if it exists.
        """
        if not (self.isparallel and self.istotal):
            raise ValueError(utils.value_err(
                self,
                'expecting a parallel pair of total maps'))

        mapping = [
                [El(n, k) for k, x in enumerate(n_data)]
                for n, n_data in enumerate(self.target.face_data)]

        to_delete = GrSet()
        shift = [[0 for _ in n_data] for n_data in mapping]

        for x in self.source:
            x1, x2 = self.fst[x], self.snd[x]
            fst, snd = mapping[x1.dim][x1.pos], mapping[x2.dim][x2.pos]
            if fst != snd:
                if fst.dim == snd.dim:
                    fst, snd = (fst, snd) if fst.pos < snd.pos else (snd, fst)
                else:
                    fst, snd = (fst, snd) if fst.dim < snd.dim else (snd, fst)
                if snd not in to_delete:
                    to_delete.add(snd)
                    for k in range(snd.pos, len(shift[snd.dim])):
                        shift[snd.dim][k] -= 1
                mapping[snd.dim][snd.pos] = mapping[fst.dim][fst.pos]
        mapping = [
                [x.shifted(shift[x.dim][x.pos]) for x in n_data]
                for n_data in mapping]

        face_data = [
                [
                    {sign:
                        {mapping[n-1][j].pos
                         for j in x[sign] if mapping[n-1][j].dim == n-1}
                     for sign in ('-', '+')}
                    for k, x in enumerate(n_data) if El(n, k) not in to_delete
                ]
                for n, n_data in enumerate(self.target.face_data)]
        quotient = OgPoset.from_face_data(face_data, wfcheck)

        return OgMap(self.target, quotient, mapping, wfcheck=False)

    def pushout(self,
                wfcheck=True):
        """
        Returns the pushout of a span of total maps, if it exists.
        """
        if not (self.isspan and self.istotal):
            raise ValueError(utils.value_err(
                self,
                'expecting a span of injective total maps'))

        coproduct = OgPoset.coproduct(self.fst.target, self.snd.target)
        coequaliser = self.then(coproduct).coequaliser(wfcheck)

        return coproduct.then(coequaliser)
