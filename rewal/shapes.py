"""
Implements shapes of cells and diagrams.
"""

from rewal import utils
from rewal.ogposets import (El, OgPoset, GrSet, GrSubset,
                            OgMap, OgMapPair)


class Shape(OgPoset):
    """
    Class for shapes of cells and diagrams.
    """
    def __new__(self):
        return OgPoset.__new__(Empty)

    # Redefining to be more lax wrt subclasses
    def __eq__(self, other):
        return isinstance(other, Shape) and \
                self.face_data == other.face_data and \
                self.coface_data == other.coface_data

    def __mul__(self, other):
        return Shape.gray(self, other)

    def __pow__(self, other):
        utils.typecheck(other, {'type': int})
        return Shape.gray(*[self for _ in range(other)])

    @property
    def isatom(self):
        """
        Returns whether the shape is an atom (has a greatest element).
        """
        return len(self.maximal()) == 1

    @property
    def isround(self):
        """
        Returns whether the shape is round (has spherical boundary).
        The result is stored after the first run.
        """
        boundary_in = self.all().boundary('-')
        boundary_out = self.all().boundary('+')
        intersection = boundary_in.intersection(boundary_out)
        for k in range(self.dim-2, -1, -1):
            boundary_in = boundary_in.boundary('-')
            boundary_out = boundary_out.boundary('+')
            if not intersection.issubset(boundary_in.union(boundary_out)):
                return False
            intersection = boundary_in.intersection(boundary_out)
        return True

    # Main constructors
    @staticmethod
    def atom(fst, snd):
        """
        Given two shapes with identical round boundaries, returns a new
        atomic shape whose input boundary is the first one and output
        boundary the second one.
        """
        for u in fst, snd:
            utils.typecheck(u, {
                'type': Shape,
                'st': lambda v: v.isround,
                'why': 'expecting a round Shape'})
        if fst.dim != snd.dim:
            raise ValueError(utils.value_err(
                snd, 'dimension does not match dimension of {}'.format(
                    repr(fst))))
        dim = fst.dim

        if dim == -1:  # Check explicitly for some simple cases.
            return Point()
        if dim == 0:
            return Arrow()
        if isinstance(fst, Globe) and isinstance(snd, Globe):
            return Shape.suspend(fst)

        in_span = OgMapPair(
                fst.boundary_inclusion('-'), snd.boundary_inclusion('-'))
        out_span = OgMapPair(
                fst.boundary_inclusion('+'), snd.boundary_inclusion('+'))
        if not in_span.isspan:
            raise ValueError(utils.value_err(
                snd, 'input boundary does not match '
                'input boundary of {}'.format(repr(fst))))
        if not out_span.isspan:
            raise ValueError(utils.value_err(
                snd, 'output boundary does not match '
                'output boundary of {}'.format(repr(fst))))

        glue_in = in_span.pushout()
        glue_out = out_span.then(glue_in).coequaliser()
        inclusion = glue_in.then(glue_out)

        sphere = inclusion.target
        # Add a greatest element
        face_data = sphere.face_data + [
                [
                    {'-':
                        {inclusion.fst[x].pos for x in fst[dim]},
                     '+':
                        {inclusion.snd[x].pos for x in snd[dim]}}
                ]]
        coface_data = sphere.coface_data + [
                [{'-': set(), '+': set()}]]
        for x in fst[dim]:
            coface_data[dim][inclusion.fst[x].pos]['-'].add(0)
        for x in snd[dim]:
            coface_data[dim][inclusion.snd[x].pos]['+'].add(0)

        new_atom = Shape.__reorder(
                OgPoset(face_data, coface_data,
                        wfcheck=False, matchcheck=False)).source

        # Recursion for positive opetopes
        if isinstance(fst, OpetopeTree) and isinstance(snd, Opetope):
            new_atom = Opetope._Shape__upgrade(new_atom)

        return new_atom

    @staticmethod
    def paste(fst, snd, dim=None):
        """
        Returns the pasting of two shapes along the output k-boundary
        of the first and the input k-boundary of the second.
        """
        if dim is None:  # default is principal composition
            dim = min(fst.dim, snd.dim) - 1
        for u in fst, snd:
            utils.typecheck(u, {'type': Shape})
        utils.typecheck(dim, {
            'type': int,
            'st': lambda n: n >= 0,
            'why': 'expecting non-negative integer'})

        span = OgMapPair(
                fst.boundary_inclusion('+', dim),
                snd.boundary_inclusion('-', dim))
        if not span.isspan:
            raise ValueError(utils.value_err(
                    snd,
                    'input {}-boundary does not match '
                    'output {}-boundary of {}'.format(
                        str(dim), str(dim), repr(fst))))

        if dim >= fst.dim:
            return snd
        if dim >= snd.dim:
            return fst
        pushout = span.pushout()
        paste = Shape.__reorder(pushout.target).source

        # Upgrade to named shape classes
        if isinstance(fst, Theta) and isinstance(snd, Theta):
            if isinstance(fst, GlobeString) and isinstance(snd, GlobeString) \
                    and fst.dim == snd.dim == dim+1:
                return GlobeString._Shape__upgrade(paste)
            return Theta._Shape__upgrade(paste)

        if isinstance(fst, OpetopeTree) and isinstance(snd, GlobeString) \
                and fst.dim == snd.dim == dim+1:
            return OpetopeTree._Shape__upgrade(paste)

        return paste

    # Other constructors
    @staticmethod
    def suspend(shape, n=1):
        utils.typecheck(shape, {'type': Shape})
        if n == 0:
            return shape

        # Suspensions of Empty are not shapes
        if isinstance(shape, Empty):
            return OgPoset.suspend(shape, n)
        if isinstance(shape, Point) and n == 1:
            return Arrow()

        suspension = Shape.__reorder(
                OgPoset.suspend(shape, n)).source

        # Theta, GlobeString, Globe are closed under suspension
        if isinstance(shape, Theta):
            if isinstance(shape, GlobeString):
                if isinstance(shape, Globe):
                    return Globe._Shape__upgrade(suspension)
                return GlobeString._Shape__upgrade(suspension)
            return Theta._Shape__upgrade(suspension)
        return suspension

    @staticmethod
    def gray(*shapes):
        """ Gray product of shapes. """
        for x in shapes:
            utils.typecheck(x, {'type': Shape})
        if len(shapes) == 0:
            return Point()
        if len(shapes) == 1:
            return shapes[0]

        oggray = OgPoset.gray(*shapes)
        if oggray in shapes:
            return oggray

        gray = Shape.__reorder(oggray).source
        # Cubes are closed under Gray products
        if all([isinstance(x, Cube) for x in shapes]):
            return Cube._Shape__upgrade(gray)
        return gray

    # Named shapes
    @staticmethod
    def empty():
        return Empty()

    @staticmethod
    def point():
        return Point()

    @staticmethod
    def arrow():
        return Arrow()

    @staticmethod
    def cube(dim=0):
        utils.typecheck(dim, {
            'type': int,
            'st': lambda n: n >= 0,
            'why': 'expecting non-negative integer'})
        arrow = Arrow()
        return Shape.gray(*[arrow for _ in range(dim)])

    @staticmethod
    def globe(dim=0):
        """ Globes. """
        return Shape.suspend(Point(), dim)

    @staticmethod
    def theta(*thetas):
        """ Batanin cells. """
        if len(thetas) > 0:
            theta = thetas[0]
            utils.typecheck(theta, {'type': Theta})
            thetas = thetas[1:]
            if len(thetas) > 0:
                return Shape.paste(
                        Shape.suspend(theta),
                        Shape.theta(*thetas), 0)
            return Shape.suspend(theta)
        return Point()

    # Special maps
    def id(self):
        return ShapeMap(super().id(),
                        wfcheck=False)

    def boundary_inclusion(self, sign=None, dim=None):
        """
        Input and output boundaries of Shapes are Shapes.
        """
        if isinstance(dim, int) and dim >= self.dim:
            return self.id()

        boundary_ogmap = super().boundary_inclusion(sign, dim)
        if sign is None:
            return boundary_ogmap
        reordering = Shape.__reorder(boundary_ogmap.source)

        return ShapeMap(reordering.then(boundary_ogmap),
                        wfcheck=False)

    def under(self, element):
        """
        Returns the inclusion of an atom as the closure of an element
        in the shape.
        """
        utils.typecheck(element, {
            'type': El,
            'st': lambda x: x in self,
            'why': 'not an element'})
        underset = GrSubset(
                GrSet(element), self,
                wfcheck=False).closure()
        reordering = Shape.__reorder(underset.as_map.source)

        return ShapeMap(reordering.then(underset.as_map),
                        wfcheck=False)

    def initial(self):
        """
        Returns the unique map from the initial (empty) shape.
        """
        return ShapeMap(
                OgMap(Empty(), self,
                      wfcheck=False),
                wfcheck=False)

    def terminal(self):
        """
        Returns the unique map to the point.
        """
        mapping = [
                [El(0, 0) for _ in n_data]
                for n_data in self.face_data]
        return ShapeMap(
                OgMap(self, Point(), mapping,
                      wfcheck=False),
                wfcheck=False)

    # Private methods
    @staticmethod
    def __reorder(shape):
        """
        Traverses all elements and returns an isomorphism
        from the shape with elements reordered in traversal order.
        """
        mapping = [[] for _ in range(shape.dim + 1)]
        marked = shape.none()

        focus_stack = [shape.all()]  # traversal begins
        while len(focus_stack) > 0:
            focus = focus_stack[-1]
            if focus.issubset(marked):
                del focus_stack[-1]
            else:
                dim = focus.dim
                if dim == 0:
                    for x in focus[dim]:
                        mapping[dim].append(x)
                        marked.support.add(x)
                    del focus_stack[-1]
                else:
                    focus_input = focus.boundary('-')
                    if not focus_input.issubset(marked):
                        focus_stack.append(focus_input)
                    else:
                        if len(focus.maximal()) == 1:
                            for x in focus[dim]:
                                mapping[dim].append(x)
                                marked.support.add(x)
                            del focus_stack[-1]
                            focus_stack.append(focus.boundary('+'))
                        else:
                            def candidates(x):
                                return [y for y in shape.cofaces(x, '-')
                                        if y in focus and y not in marked]
                            x = next(
                                    x for x in mapping[dim-1]
                                    if x in focus and len(candidates(x)) > 0)
                            focus_stack.append(GrSubset(
                                GrSet(candidates(x)[0]),
                                shape, wfcheck=False).closure())

        def reordered_faces(x, sign):
            return {k for y in shape.faces(x, sign)
                    for k in range(shape.size[x.dim - 1])
                    if y == mapping[x.dim - 1][k]}

        def reordered_cofaces(x, sign):
            return {k for y in shape.cofaces(x, sign)
                    for k in range(shape.size[x.dim + 1])
                    if y == mapping[x.dim + 1][k]}

        face_data = [
                [
                    {sign: reordered_faces(x, sign)
                     for sign in ('-', '+')}
                    for x in n_data
                ]
                for n_data in mapping]
        coface_data = [
                [
                    {sign: reordered_cofaces(x, sign)
                     for sign in ('-', '+')}
                    for x in n_data
                ]
                for n_data in mapping]
        reordered = Shape.__upgrade(
                OgPoset(face_data, coface_data,
                        wfcheck=False, matchcheck=False))

        return OgMap(reordered, shape, mapping,
                     wfcheck=False)

    @classmethod
    def __upgrade(cls, ogp):
        """
        Forces upgrade of an OgPoset to the shape class.
        """
        shape = OgPoset.__new__(cls)
        shape._face_data = ogp.face_data
        shape._coface_data = ogp.coface_data
        return shape

    @classmethod
    def __upgrademapsrc(cls, ogmap):
        """
        Upgrades the source of a map to the shape class.
        """
        return OgMap(
                cls._Shape__upgrade(ogmap.source),
                ogmap.target,
                ogmap.mapping,
                wfcheck=False)

    @classmethod
    def __upgrademaptgt(cls, ogmap):
        """
        Upgrades the target of a map to the shape class.
        """
        return OgMap(
                ogmap.source,
                cls._Shape__upgrade(ogmap.target),
                ogmap.mapping,
                wfcheck=False)


class Simplex(Shape):
    """
    Class for oriented simplices.
    """
    def __new__(self):
        return OgPoset.__new__(Empty)


class Empty(Simplex):
    """
    Class for the empty shape.
    """
    def __new__(self):
        return OgPoset.__new__(Empty)

    def __init__(self):
        super().__init__([], [],
                         wfcheck=False, matchcheck=False)


class Cube(Shape):
    """
    Class for cubes.
    """
    def __new__(self):
        return OgPoset.__new__(Point)

    def cube_face(self, k, sign):
        """ Cubical face maps. """
        sign = utils.mksign(sign)
        utils.typecheck(k, {
            'type': int,
            'st': lambda k: k in range(self.dim),
            'why': 'out of bounds'})
        basic_faces = {
            '-': ShapeMap(OgMap(
                Point(), Arrow(),
                [[El(0, 0)]],
                wfcheck=False), wfcheck=False),
            '+': ShapeMap(OgMap(
                Point(), Arrow(),
                [[El(0, 1)]],
                wfcheck=False), wfcheck=False)}
        arrowid = Arrow().id()
        maps = [arrowid for _ in range(k)] + \
            [basic_faces[sign]] + \
            [arrowid for _ in range(self.dim - k - 1)]
        return ShapeMap.gray(*maps)

    def cube_degeneracy(self, k):
        """ Cubical degeneracy maps. """
        utils.typecheck(k, {
            'type': int,
            'st': lambda k: k in range(self.dim),
            'why': 'out of bounds'})
        arrowid = Arrow().id()
        maps = [arrowid for _ in range(k)] + \
            [Arrow().terminal()] + \
            [arrowid for _ in range(self.dim - k - 1)]
        return ShapeMap.gray(*maps)

    def connection(self, k, sign):
        """ Cubical connection maps. """
        sign = utils.mksign(sign)
        utils.typecheck(k, {
            'type': int,
            'st': lambda k: k in range(self.dim),
            'why': 'out of bounds'})
        basic_connections = {
                '-': ShapeMap(OgMap(
                        Shape.cube(2), Arrow(),
                        [[El(0, 0), El(0, 0), El(0, 1), El(0, 0)],
                         [El(0, 0), El(1, 0), El(0, 0), El(1, 0)],
                         [El(1, 0)]],
                        wfcheck=False), wfcheck=False),
                '+': ShapeMap(OgMap(
                        Shape.cube(2), Arrow(),
                        [[El(0, 0), El(0, 1), El(0, 1), El(0, 1)],
                         [El(1, 0), El(0, 1), El(1, 0), El(0, 1)],
                         [El(1, 0)]],
                        wfcheck=False), wfcheck=False)}
        arrowid = Arrow().id()
        maps = [arrowid for _ in range(k)] + \
            [basic_connections[sign]] + \
            [arrowid for _ in range(self.dim - k - 1)]
        return ShapeMap.gray(*maps)


class Theta(Shape):
    """
    Class for the objects of the Theta category (Batanin cells).
    """
    def __new__(self):
        return OgPoset.__new__(Point)


class OpetopeTree(Shape):
    """
    Class for shapes that can appear as input boundaries of positive
    opetopes.
    """
    def __new__(self):
        return OgPoset.__new__(Point)


class GlobeString(Theta, OpetopeTree):
    """
    Class for 'strings of globes' in the top dimension, the intersection
    of opetope trees and Batanin cells.
    """
    def __new__(self):
        return OgPoset.__new__(Point)


class Opetope(OpetopeTree):
    """
    Class for positive opetopes.
    """
    def __new__(self):
        return OgPoset.__new__(Point)


class Globe(GlobeString, Opetope):
    """ Class for globes. """
    def __new__(self):
        return OgPoset.__new__(Point)


class Point(Globe, Simplex, Cube):
    """ Class for the point. """
    def __new__(self):
        return OgPoset.__new__(Point)

    def __init__(self):
        super().__init__(
                [[{'-': set(), '+': set()}]],
                [[{'-': set(), '+': set()}]],
                wfcheck=False, matchcheck=False)


class Arrow(Globe, Simplex, Cube):
    """ Class for the arrow. """
    def __new__(self):
        return OgPoset.__new__(Arrow)

    def __init__(self):
        super().__init__(
                [
                    [{'-': set(), '+': set()}, {'-': set(), '+': set()}],
                    [{'-': {0}, '+': {1}}]
                ], [
                    [{'-': {0}, '+': set()}, {'-': set(), '+': {0}}],
                    [{'-': set(), '+': set()}]
                ],
                wfcheck=False, matchcheck=False)


class ShapeMap(OgMap):
    """
    An OgMap overlay for total maps between Shape objects.
    Used for constructions that are not well-defined for general
    maps between oriented graded posets.
    """

    def __init__(self, ogmap, wfcheck=True):
        if wfcheck:
            utils.typecheck(ogmap, {'type': OgMap})
            for x in ogmap.source, ogmap.target:
                utils.typecheck(x, {'type': Shape})
            if not ogmap.istotal:
                raise ValueError(utils.value_err(
                    ogmap,
                    'a ShapeMap must be total'))

        super().__init__(ogmap.source, ogmap.target, ogmap.mapping,
                         wfcheck=False)

    def __mul__(self, other):
        return ShapeMap.gray(self, other)

    def __pow__(self, other):
        utils.typecheck(other, {'type': int})
        return ShapeMap.gray(*[self for _ in range(other)])

    @staticmethod
    def gray(*maps):
        for f in maps:
            utils.typecheck(f, {'type': ShapeMap})
        if len(maps) == 0:
            return Point().id()
        if len(maps) == 1:
            return maps[0]

        def oggray(fst, snd, *others):
            if len(others) > 0:
                return oggray(oggray(fst, snd), *others)

            size1 = fst.target.size + [0 for _ in range(snd.target.dim)]
            size2 = snd.target.size + [0 for _ in range(fst.target.dim)]

            def pair(x, y):
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

        gray = oggray(*maps)
        if gray.source in [f.source for f in maps]:
            if gray.target in [f.target for f in maps]:
                return gray
        if gray.source not in [f.source for f in maps]:
            gray = Shape._Shape__reorder(gray.source).then(gray)
            if all([isinstance(x, Cube) for x in [f.source for f in maps]]):
                gray = Cube._Shape__upgrademapsrc(gray)
        if gray.target not in [f.target for f in maps]:
            gray = gray.then(Shape._Shape__reorder(gray.target).inv())
            if all([isinstance(x, Cube) for x in [f.target for f in maps]]):
                gray = Cube._Shape__upgrademaptgt(gray)

        return ShapeMap(gray,
                        wfcheck=False)


def atom(fst, snd):
    return Shape.atom(fst, snd)


def paste(fst, snd, dim=None):
    return Shape.paste(fst, snd, dim)


def suspend(shape, n=1):
    return Shape.suspend(shape, n)
