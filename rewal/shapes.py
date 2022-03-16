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

    _isround = None

    def __init__(self):
        """
        Initialises to the empty shape.
        """
        self._face_data = []
        self._coface_data = []
        self._isatom = False

    @property
    def isatom(self):
        """
        Returns whether the shape is an atom (has a greatest element).
        We let this be stored at construction.
        """
        return self._isatom

    @property
    def isround(self):
        """
        Returns whether the shape has spherical boundary. The result
        is stored after the first run.
        """
        if self._isround is not None:
            return self._isround
        boundary_in = self.all().boundary('-')
        boundary_out = self.all().boundary('+')
        intersection = boundary_in.intersection(boundary_out)
        for k in range(self.dim-2, -1, -1):
            boundary_in = boundary_in.boundary('-')
            boundary_out = boundary_out.boundary('+')
            if not intersection.issubset(boundary_in.union(boundary_out)):
                self._isround = False
                return False
            intersection = boundary_in.intersection(boundary_out)
        self._isround = True
        return True

    # Constructors
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
                'why': 'expecting a shape with spherical boundary'})
        if fst.dim != snd.dim:
            raise ValueError(utils.value_err(
                snd, 'dimension does not match dimension of {}'.format(
                    repr(fst))))
        dim = fst.dim

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

        new_atom = inclusion.target
        # Add a greatest element
        new_atom.face_data.append(
                [
                    {'-':
                        {inclusion.fst[x].pos for x in fst[dim]},
                     '+':
                        {inclusion.snd[x].pos for x in snd[dim]}}
                ])
        new_atom.coface_data.append([{'-': set(), '+': set()}])
        for x in fst[dim]:
            new_atom.coface_data[dim][inclusion.fst[x].pos]['-'].add(0)
        for x in snd[dim]:
            new_atom.coface_data[dim][inclusion.snd[x].pos]['+'].add(0)

        return Shape.__reorder(new_atom).source

    @staticmethod
    def paste_cospan(fst, snd, dim=None):
        """
        Returns the pasting of two shapes along their dim-boundary,
        together with their inclusions into the pasting.
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
            return OgMapPair(fst.id(), span.fst)
        if dim >= snd.dim:
            return OgMapPair(span.snd, snd.id())
        pushout = span.pushout()
        reordering = Shape.__reorder(pushout.target).inv()
        paste_cospan = pushout.then(reordering)

        return OgMapPair(
                ShapeMap(paste_cospan.fst, wfcheck=False),
                ShapeMap(paste_cospan.snd, wfcheck=False))

    @staticmethod
    def paste(fst, snd, dim=None):
        return Shape.paste_cospan(fst, snd, dim).target

    # Some special maps of shapes.
    def id(self):
        return ShapeMap(super().id(),
                        wfcheck=False)

    def boundary_inclusion(self, sign=None, dim=None):
        """
        Input and output boundaries of Shapes are Shapes.
        """
        boundary_ogmap = super().boundary_inclusion(sign, dim)
        if sign is None:
            return boundary_ogmap
        reordering = Shape.__reorder(boundary_ogmap.source)
        return ShapeMap(reordering.then(boundary_ogmap))

    def initial(self):
        """
        Returns the unique map from the initial (empty) shape.
        """
        return ShapeMap(
                OgMap(Shape(), self, wfcheck=False),
                wfcheck=False)

    def terminal(self):
        """
        Returns the unique map to the terminal shape (the point).
        """
        mapping = [
                [El(0, 0) for _ in n_data]
                for n_data in self.face_data]
        return ShapeMap(
                OgMap(self, Shape.point(), mapping, wfcheck=False),
                wfcheck=False)

    # Some special named shapes.
    @staticmethod
    def globe(dim):  # TODO: make this use suspensions
        """ The globes. """
        utils.typecheck(dim, {'type': int})
        if dim >= 0:
            lower = Shape.globe(dim - 1)
            return Shape.atom(lower, lower)
        return Shape()

    @classmethod
    def point(cls):
        """ The point. """
        face_data = [[{'-': set(), '+': set()}]]
        coface_data = [[{'-': set(), '+': set()}]]
        return cls.__upgrade(
                OgPoset(face_data, coface_data, wfcheck=False))

    @classmethod
    def arrow(cls):
        """ The arrow. """
        return cls.__upgrade(Shape.globe(1))

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
                focus_input = focus.boundary('-')
                if not focus_input.issubset(marked):
                    focus_stack.append(focus_input)
                else:
                    if len(focus.maximal()) == 1:
                        for x in focus.maximal():
                            mapping[x.dim].append(x)
                            marked.support.add(x)
                        del focus_stack[-1]
                        focus_stack.append(focus.boundary('+'))
                    else:
                        def candidates(x):
                            return [y for y in shape.cofaces(x, '-')
                                    if y not in marked]
                        x = next(
                                x for x in mapping[focus.support.dim - 1]
                                if len(candidates(x)) > 0)
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
                        wfcheck=False))

        return OgMap(reordered, shape, mapping,
                     wfcheck=False)

    @classmethod
    def __upgrade(cls, ogposet):
        """
        Forces upgrade of an OgPoset to the shape class.
        """
        shape = cls.__new__(cls)
        shape._face_data = ogposet.face_data
        shape._coface_data = ogposet.coface_data
        shape._isatom = True if len(ogposet.maximal()) == 1 else False
        return shape


class ShapeMap(OgMap):
    """
    An OgMap overlay for total maps between Shape objects.
    Used for constructions that are not well-defined for general
    maps between oriented graded posets.
    """

    _isinjective = None
    _issurjective = None

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

    @property
    def isinjective(self):
        """ We will store the result after the first run. """
        if self._isinjective is None:
            self._isinjective = super().isinjective
        return self._isinjective

    @property
    def issurjective(self):
        """ We will store the result after the first run. """
        if self._issurjective is None:
            self._issurjective = super().issurjective
        return self._issurjective


def atom(fst, snd):
    return Shape.atom(fst, snd)


def paste(fst, snd, dim=None):
    return Shape.paste(fst, snd, dim)


def globe(dim):
    return Shape.globe(dim)
