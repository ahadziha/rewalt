"""
Implements shapes of cells and diagrams.
"""

from rewal import utils
from rewal.ogposets import (El, OgPoset, GrSet, GrSubset, Closed,
                            OgMap, OgMapPair)


class Shape(OgPoset):
    """
    Class for shapes of cells and diagrams.
    """
    def __new__(self):
        return OgPoset.__new__(Empty)

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
        """
        return self.all().isround

    # Main constructors
    @staticmethod
    def atom_cospan(fst, snd):
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

        if dim == -1:  # Avoid more work in this simple case
            return OgMapPair(
                    Shape.point().initial(),
                    Shape.point().initial()
                    )

        in_span = OgMapPair(
                fst.boundary('-'), snd.boundary('-'))
        out_span = OgMapPair(
                fst.boundary('+'), snd.boundary('+'))
        if not in_span.isspan:
            raise ValueError(utils.value_err(
                snd, 'input boundary does not match '
                'input boundary of {}'.format(repr(fst))))
        if not out_span.isspan:
            raise ValueError(utils.value_err(
                snd, 'output boundary does not match '
                'output boundary of {}'.format(repr(fst))))

        glue_in = in_span.pushout(wfcheck=False)
        glue_out = out_span.then(glue_in).coequaliser(wfcheck=False)
        inclusion = glue_in.then(glue_out)

        sphere = inclusion.target
        # Add a greatest element
        face_data = [
                *sphere.face_data,
                [
                    {'-':
                        {inclusion.fst[x].pos for x in fst[dim]},
                     '+':
                        {inclusion.snd[x].pos for x in snd[dim]}}
                ]]
        coface_data = [
                *sphere.coface_data,
                [{'-': set(), '+': set()}]]
        for x in fst[dim]:
            coface_data[dim][inclusion.fst[x].pos]['-'].add(0)
        for x in snd[dim]:
            coface_data[dim][inclusion.snd[x].pos]['+'].add(0)

        new_atom = OgPoset(face_data, coface_data,
                           wfcheck=False, matchcheck=False)
        boundary_in = OgMap(
            fst, new_atom, inclusion.fst.mapping,
            wfcheck=False)
        boundary_out = OgMap(
            snd, new_atom, inclusion.snd.mapping,
            wfcheck=False)

        atom_cospan = OgMapPair(boundary_in, boundary_out).then(
                Shape._reorder(new_atom).inv())

        def inheritance():
            if isinstance(fst, OpetopeTree) and isinstance(snd, Opetope):
                if isinstance(fst, Globe):
                    if fst.dim == 0:
                        return Arrow
                    return Globe
                return Opetope
            return Shape

        return inheritance()._upgrademaptgt(atom_cospan)

    @staticmethod
    def atom(fst, snd):
        return Shape.atom_cospan(fst, snd).target

    def _atom(self, other):
        return Shape.atom(self, other)

    @staticmethod
    def paste_cospan(fst, snd, dim=None):
        """
        Returns the pasting of two shapes along the output k-boundary
        of the first and the input k-boundary of the second.
        """
        for u in fst, snd:
            utils.typecheck(u, {'type': Shape})
        if dim is None:  # default is principal composition
            dim = min(fst.dim, snd.dim) - 1
        utils.typecheck(dim, {
            'type': int,
            'st': lambda n: n >= 0,
            'why': 'expecting non-negative integer'})

        span = OgMapPair(
                fst.boundary('+', dim),
                snd.boundary('-', dim))
        if not span.isspan:
            raise ValueError(utils.value_err(
                    snd,
                    'input {}-boundary does not match '
                    'output {}-boundary of {}'.format(
                        str(dim), str(dim), repr(fst))))

        if dim >= fst.dim:
            return OgMapPair(
                    span.snd, snd.id())
        if dim >= snd.dim:
            return OgMapPair(
                    fst.id(), span.fst)

        pushout = span.pushout(wfcheck=False)
        pasting = pushout.then(Shape._reorder(pushout.target).inv())

        def inheritance():
            if isinstance(fst, Theta) and isinstance(snd, Theta):
                if isinstance(fst, GlobeString) and \
                        isinstance(snd, GlobeString) \
                        and fst.dim == snd.dim == dim+1:
                    return GlobeString
                return Theta
            if isinstance(fst, OpetopeTree) and isinstance(snd, GlobeString) \
                    and fst.dim == snd.dim == dim+1:
                return OpetopeTree
            return Shape

        return inheritance()._upgrademaptgt(pasting)

    @staticmethod
    def paste(fst, snd, dim=None):
        return Shape.paste_cospan(fst, snd, dim).target

    def _paste(self, other, dim=None):
        return Shape.paste(self, other, dim)

    # Other constructors
    @staticmethod
    def paste_along_cospan(fst, snd):
        """
        Returns the pasting of two shapes along the entire input
        (output) k-boundary of one, and a subshape of the output
        (input) k-boundary of the other.
        """
        span = OgMapPair(fst, snd)
        utils.typecheck(span, {
            'type': OgMapPair,
            'st': lambda x: x.isspan and x.isinjective,
            'why': 'expecting a span of injective maps'
            }, {'type': ShapeMap})

        dim = span.source.dim
        fst_image = fst.image()
        snd_image = snd.image()
        fst_output = fst.target.all().boundary('+', dim)
        snd_input = snd.target.all().boundary('-', dim)

        def condition():
            t1 = fst_image.issubset(fst_output)
            t2 = snd_image.issubset(snd_input)
            t3 = fst_image == fst_output or snd_image == snd_input
            return t1 and t2 and t3

        if not condition():
            raise ValueError(utils.value_err(
                span, 'not a well-formed span for pasting'))
        if fst_image == fst_output:
            if snd_image == snd_input:
                return Shape.paste_cospan(fst.target, snd.target, dim)
            if not Shape._issubmol(snd_image, snd_input):
                raise ValueError(utils.value_err(
                    snd, 'cannot paste along this map'))
        else:
            if not Shape._issubmol(fst_image, fst_output):
                raise ValueError(utils.value_err(
                    fst, 'cannot paste along this map'))

        pushout = span.pushout(wfcheck=False)
        pasting = pushout.then(Shape._reorder(pushout.target).inv())

        def inheritance():
            if isinstance(fst.target, OpetopeTree) and \
                    isinstance(snd.target, OpetopeTree) \
                    and fst.target.dim == snd.target.dim == dim+1:
                return OpetopeTree
            return Shape

        return inheritance()._upgrademaptgt(pasting)

    @staticmethod
    def paste_along(fst, snd):
        return Shape.paste_along_cospan(fst, snd).target

    def paste_at_output_cospan(self, pos, other):
        """
        Paste along the inclusion of one of the outputs.
        """
        return Shape.paste_along_cospan(
            self.atom_inclusion(El(self.dim-1, pos)),
            other.boundary('-'))

    def paste_at_output(self, pos, other):
        return self.paste_at_output_cospan(pos, other).target

    def paste_at_input_cospan(self, pos, other):
        """
        Paste along the inclusion of one of the inputs.
        """
        return Shape.paste_along(
            other.boundary('+'),
            self.atom_inclusion(El(self.dim-1, pos)))

    def paste_at_input(self, pos, other):
        return self.paste_at_input_cospan(pos, other).target

    @staticmethod
    def suspend(shape, n=1):
        if n == 0:
            return shape

        if not isinstance(shape, Shape) or isinstance(shape, Empty):
            return OgPoset.suspend(shape, n)
        if isinstance(shape, Point) and n == 1:
            return Arrow()

        suspension = Shape._reorder(
                OgPoset.suspend(shape, n)).source

        def inheritance():
            if isinstance(shape, Theta):
                if isinstance(shape, GlobeString):
                    if isinstance(shape, Globe):
                        return Globe
                    return GlobeString
                return Theta
            return Shape

        return inheritance()._upgrade(suspension)

    @staticmethod
    def gray(*shapes):
        """ Gray product of shapes. """
        for x in shapes:
            if not isinstance(x, Shape):
                return OgPoset.gray(*shapes)
        if len(shapes) == 0:
            return Point()
        if len(shapes) == 1:
            return shapes[0]

        oggray = OgPoset.gray(*shapes)
        if oggray in shapes:
            return oggray

        gray = Shape._reorder(oggray).source

        def inheritance(l):
            if all([isinstance(x, Cube) for x in l]):
                return Cube
            return Shape

        return inheritance(shapes)._upgrade(gray)

    @staticmethod
    def join(*shapes):
        for x in shapes:
            if not isinstance(x, Shape):
                return OgPoset.join(*shapes)
        if len(shapes) == 0:
            return Empty()
        if len(shapes) == 1:
            return shapes[0]

        ogjoin = OgPoset.join(*shapes)
        if len(ogjoin) == 3:  # check explicitly if it's arrow
            return Arrow()
        if ogjoin in shapes:
            return ogjoin

        join = Shape._reorder(ogjoin).source

        def inheritance(l):
            if all([isinstance(x, Simplex) for x in l]):
                return Simplex
            return Shape

        return inheritance(shapes)._upgrade(join)

    @staticmethod
    def dual(shape, *dims):
        dual = Shape._reorder(OgPoset.dual(shape, *dims)).source
        if shape == dual:
            return shape

        def inheritance():
            if isinstance(shape, Theta):
                return Theta
            return Shape

        return inheritance()._upgrade(dual)

    def merge(self):
        """
        Returns an atom with the same boundary as the shape (if round).
        """
        if self.isatom:
            return self
        if not self.isround:
            raise ValueError(utils.value_err(
                self, 'not a round shape'))

        merged = Shape.atom(
                self.boundary('-').source,
                self.boundary('+').source)

        def inheritance():
            if isinstance(self, OpetopeTree):
                if isinstance(self, GlobeString):
                    if self.dim == 1:
                        return Arrow
                    return Globe
                return Opetope
            return Shape

        return inheritance()._upgrade(merged)

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
    def simplex(dim=-1):
        utils.typecheck(dim, {
            'type': int,
            'st': lambda n: n >= -1,
            'why': 'expecting integer >= -1'})
        point = Point()
        return Shape.join(*[point for _ in range(dim+1)])

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

    def boundary(self, sign=None, dim=None):
        """
        Input and output boundaries of Shapes are Shapes.
        """
        if isinstance(dim, int) and dim >= self.dim:
            return self.id()

        boundary_ogmap = super().boundary(sign, dim)
        if sign is None:
            return boundary_ogmap

        reordering = Shape._reorder(boundary_ogmap.source)
        boundary = reordering.then(boundary_ogmap)

        def inheritance():
            if isinstance(self, OpetopeTree):
                if isinstance(self, GlobeString):
                    return Globe
                if utils.mksign(sign) == '+':
                    return Opetope
                return OpetopeTree
            return Shape

        return inheritance()._upgrademapsrc(boundary)

    def atom_inclusion(self, element):
        """
        Returns the inclusion of an atom as the closure of an element
        in the shape.
        """
        utils.typecheck(element, {
            'type': El,
            'st': lambda x: x in self,
            'why': 'not an element of {}'.format(repr(self))})
        underset = GrSubset(
                GrSet(element), self,
                wfcheck=False).closure()
        reordering = Shape._reorder(underset.as_map.source)
        inclusion = reordering.then(underset.as_map)

        def inheritance():
            if element.dim == 0:
                return Point
            if element.dim == 1:
                return Arrow
            if isinstance(self, Theta):
                return Globe
            if isinstance(self, OpetopeTree):
                return Opetope
            if isinstance(self, Simplex):
                return Simplex
            if isinstance(self, Cube):
                return Cube

        return inheritance()._upgrademapsrc(inclusion)

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

    def inflate(self, collapsed=None):
        """
        Forms a cylinder on a shape with some "collapsed" sides, specified
        by a closed subset of the boundary of the shape, and returns its
        projection on the original shape.
        Used in constructing units and unitors on diagrams.
        """
        if self.dim == -1:  # Some simple cases
            return self.id()
        if self.dim == 0:
            return Shape.arrow().terminal()

        boundary_set = self.all().boundary()
        if collapsed is not None:
            utils.typecheck(collapsed, {
                'type': Closed,
                'st': lambda x: x.issubset(boundary_set),
                'why': "expecting a closed subset of the shape's boundary"})
        else:
            collapsed = boundary_set  # Default is whole boundary.

        asmap = collapsed.as_map
        arrow = Shape.arrow()
        map1 = OgMap.gray(arrow.id(), asmap)
        map2 = OgMap.gray(arrow.terminal(), asmap.source.id())
        pushout = OgMapPair(map1, map2).pushout(wfcheck=False)

        # We use the cylinder projection map and the first leg of the pushout
        # to define the projection map.
        cyl_proj = OgMap.gray(arrow.terminal(), self.id())
        collapse = pushout.fst
        mapping = [
                [None for _ in n_data]
                for n_data in pushout.target.face_data]
        for x in collapse.source:
            y = collapse[x]
            mapping[y.dim][y.pos] = cyl_proj[x]

        ogproj = OgMap(collapse.target, self, mapping,
                       wfcheck=False)
        proj = Shape._reorder(collapse.target).then(ogproj)

        def inheritance():
            if collapsed == boundary_set and isinstance(self, Opetope):
                if isinstance(self, Globe):
                    return Globe
                return Opetope
            return Shape

        return inheritance()._upgrademapsrc(proj)

    # Private methods
    @staticmethod
    def _reorder(shape):
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
        reordered = Shape._upgrade(
                OgPoset(face_data, coface_data,
                        wfcheck=False, matchcheck=False))

        return OgMap(reordered, shape, mapping,
                     wfcheck=False)

    @staticmethod
    def _issubmol(fst, snd):
        """
        Given two molecules in a shape (as closed subsets), returns
        whether the first one is a round submolecule of the second.
        """
        if not fst.isround:
            return False
        if fst == snd:
            return True
        if len(fst.maximal()) == 1:
            return True
        raise NotImplementedError(
            'Can only paste along atoms or the entire boundary for now.')

    @classmethod
    def _upgrade(cls, ogp):
        """
        Forces upgrade of an OgPoset to the shape class.
        """
        shape = OgPoset.__new__(cls)
        shape.__init__(ogp.face_data, ogp.coface_data,
                       wfcheck=False, matchcheck=False)
        return shape

    @classmethod
    def _upgrademapsrc(cls, ogmap):
        """
        Upgrades the source of an OgMap to a shape class, and declares
        it a ShapeMap.
        """
        if isinstance(ogmap, OgMapPair):
            return OgMapPair(
                    cls._upgrademapsrc(ogmap.fst),
                    cls._upgrademapsrc(ogmap.snd))
        return ShapeMap(OgMap(
                cls._upgrade(ogmap.source),
                ogmap.target,
                ogmap.mapping,
                wfcheck=False), wfcheck=False)

    @classmethod
    def _upgrademaptgt(cls, ogmap):
        """
        Upgrades the target of an OgMap to a shape class, and declares
        it a ShapeMap.
        """
        if isinstance(ogmap, OgMapPair):
            return OgMapPair(
                    cls._upgrademaptgt(ogmap.fst),
                    cls._upgrademaptgt(ogmap.snd))
        return ShapeMap(OgMap(
            ogmap.source,
            cls._upgrade(ogmap.target),
            ogmap.mapping,
            wfcheck=False), wfcheck=False)


class Simplex(Shape):
    """
    Class for oriented simplices.
    """
    def __new__(self):
        return OgPoset.__new__(Empty)

    def simplex_face(self, k):
        """ Simplicial face maps. """
        utils.typecheck(k, {
            'type': int,
            'st': lambda k: k in range(self.dim + 1),
            'why': 'out of bounds'})
        pointid = Point().id()
        maps = [
                *[pointid for _ in range(k)],
                Empty().terminal(),
                *[pointid for _ in range(self.dim - k)]
            ]
        return ShapeMap.join(*maps)

    def simplex_degeneracy(self, k):
        """ Simplicial degeneracy maps. """
        utils.typecheck(k, {
            'type': int,
            'st': lambda k: k in range(self.dim + 1),
            'why': 'out of bounds'})
        pointid = Point().id()
        maps = [
                *[pointid for _ in range(k)],
                Arrow().terminal(),
                *[pointid for _ in range(self.dim - k)]
            ]
        return ShapeMap.join(*maps)


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
        maps = [
                *[arrowid for _ in range(k)],
                basic_faces[sign],
                *[arrowid for _ in range(self.dim - k - 1)]
            ]
        return ShapeMap.gray(*maps)

    def cube_degeneracy(self, k):
        """ Cubical degeneracy maps. """
        utils.typecheck(k, {
            'type': int,
            'st': lambda k: k in range(self.dim + 1),
            'why': 'out of bounds'})
        arrowid = Arrow().id()
        maps = [
                *[arrowid for _ in range(k)],
                Arrow().terminal(),
                *[arrowid for _ in range(self.dim - k)]
            ]
        return ShapeMap.gray(*maps)

    def cube_connection(self, k, sign):
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
        maps = [
                *[arrowid for _ in range(k)],
                basic_connections[sign],
                *[arrowid for _ in range(self.dim - k - 1)]
            ]
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
    Class for 'strings of globes' pasted in the top dimension, the
    intersection of opetope trees and Batanin cells.
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

    def __init__(self, face_data=None, coface_data=None,
                 wfcheck=False, matchcheck=False):
        super().__init__(
                [[{'-': set(), '+': set()}]],
                [[{'-': set(), '+': set()}]],
                wfcheck=False, matchcheck=False)


class Arrow(Globe, Simplex, Cube):
    """ Class for the arrow. """
    def __new__(self):
        return OgPoset.__new__(Arrow)

    def __init__(self, face_data=None, coface_data=None,
                 wfcheck=False, matchcheck=False):
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
    """

    def __init__(self, ogmap, wfcheck=True):
        if wfcheck:
            utils.typecheck(ogmap, {
                'type': OgMap,
                'st': lambda f: f.istotal,
                'why': 'a ShapeMap must be total'})
            for x in ogmap.source, ogmap.target:
                utils.typecheck(x, {'type': Shape})

        super().__init__(ogmap.source, ogmap.target, ogmap.mapping,
                         wfcheck=False)

    def then(self, other, *others):
        for f in [other, *others]:
            if not isinstance(f, ShapeMap):
                return super().then(other, *others)
        return ShapeMap(
                super().then(other, *others),
                wfcheck=False)

    @staticmethod
    def gray(*maps):
        for f in maps:
            if not isinstance(f, ShapeMap):
                return OgMap.gray(*maps)
        if len(maps) == 0:
            return Point().id()
        if len(maps) == 1:
            return maps[0]

        gray = OgMap.gray(*maps)
        if gray.source in [f.source for f in maps]:
            if gray.target in [f.target for f in maps]:
                return gray

        if gray.source not in [f.source for f in maps]:
            gray = Shape._reorder(gray.source).then(gray)
        if gray.target not in [f.target for f in maps]:
            gray = gray.then(Shape._reorder(gray.target).inv())

        def inheritance(l):
            if all([isinstance(x, Cube) for x in l]):
                return Cube
            return Shape

        return inheritance([f.source for f in maps])._upgrademapsrc(
                inheritance([f.target for f in maps])._upgrademaptgt(
                    gray))

    @staticmethod
    def join(*maps):
        for f in maps:
            if not isinstance(f, ShapeMap):
                return OgMap.join(*maps)
        if len(maps) == 0:
            return Empty().id()
        if len(maps) == 1:
            return maps[0]

        join = OgMap.join(*maps)
        if join.source in [f.source for f in maps]:
            if join.target in [f.target for f in maps]:
                return join

        if join.source not in [f.source for f in maps]:
            join = Shape._reorder(join.source).then(join)
        if join.target not in [f.target for f in maps]:
            join = join.then(Shape._reorder(join.target).inv())

        def inheritance(l):
            if all([isinstance(x, Simplex) for x in l]):
                if sum([x.dim+1 for x in l]) == 2:
                    return Arrow
                return Simplex
            return Shape

        return inheritance([f.source for f in maps])._upgrademapsrc(
                inheritance([f.target for f in maps])._upgrademaptgt(
                    join))

    def dual(shapemap, *dims):
        if not isinstance(shapemap, ShapeMap):
            return OgMap.dual(shapemap, *dims)

        ogdual = OgMap.dual(shapemap, *dims)
        dual = Shape._reorder(ogdual.source).then(ogdual).then(
                Shape._reorder(ogdual.target).inv())

        def inheritance(x, y):
            if x == y:
                return x.__class__
            else:
                if isinstance(x, Theta):
                    return Theta
            return Shape

        return inheritance(shapemap.source, dual.source)._upgrademapsrc(
                inheritance(shapemap.target, dual.target)._upgrademaptgt(
                    dual))
