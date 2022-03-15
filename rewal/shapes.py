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
        We let this be recorded by the constructors to avoid computing
        the maximal elements at every call.
        """
        return self._isatom

    @property
    def isround(self):
        """ Returns whether the shape has spherical boundary. """
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

    # Constructors
    @staticmethod
    def atom(fst, snd):
        for u in fst, snd:
            utils.typecheck(u, {
                'type': Shape,
                'st': lambda v: v.isround,
                'why': 'expecting a shape with spherical boundary'})
        if fst.dim != snd.dim:
            raise ValueError(utils.value_err(
                snd, 'dimension does not match dimension of {}'.format(
                    repr(fst))))

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
                        {inclusion.fst[x].pos for x in fst[fst.dim]},
                     '+':
                        {inclusion.snd[x].pos for x in snd[snd.dim]}}
                ])
        new_atom.coface_data.append([{'-': set(), '+': set()}])
        for x in fst[fst.dim]:
            new_atom.coface_data[x.dim][inclusion.fst[x].pos]['-'].add(0)
        for x in snd[snd.dim]:
            new_atom.coface_data[x.dim][inclusion.snd[x].pos]['+'].add(0)

        return Shape.__reorder(new_atom).source

    @staticmethod
    def paste(fst, snd, dim):
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
            return fst
        if dim >= snd.dim:
            return snd

        return Shape.__reorder(span.pushout().target).source

    # Some special maps of shapes.
    def id(self):
        return ShapeMap.from_ogmap(super().id())

    def boundary_inclusion(self, sign=None, dim=None):
        """
        Input and output boundaries of Shapes are Shapes.
        """
        boundary_ogmap = super().boundary_inclusion(sign, dim)
        if sign is None:
            return boundary_ogmap
        reordering = Shape.__reorder(boundary_ogmap.source)
        return ShapeMap.from_ogmap(
                reordering.then(boundary_ogmap))

    def initial(self):
        """
        Returns the unique map from the initial (empty) shape.
        """
        return ShapeMap(Shape(), self,
                        wfcheck=False)

    def terminal(self):
        """
        Returns the unique map to the terminal shape (the point).
        """
        mapping = [
                [El(0, 0) for _ in n_data]
                for n_data in self.face_data]
        return ShapeMap(self, Shape.point(), mapping,
                        wfcheck=False)

    # Some special named shapes.
    @staticmethod
    def globe(dim):
        """ The globes. """
        utils.typecheck(dim, {'type': int})
        if dim >= 0:
            lower = Shape.globe(dim - 1)
            return Shape.atom(lower, lower)
        return Shape()

    @classmethod
    def point(cls):
        """ The point. """
        return cls.__upgrade(Shape.globe(0))

    @classmethod
    def arrow(cls):
        """ The arrow. """
        return cls.__upgrade(Shape.globe(1))

    # Private methods
    @staticmethod
    def __reorder(shape):
        """
        Traverses all elements of the shape and returns an isomorphism
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
        reordered_shape = Shape.__upgrade(
                OgPoset(face_data, coface_data,
                        wfcheck=False))

        return OgMap(reordered_shape, shape, mapping,
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
    An OgMap overlay for maps between Shape objects.
    Used for constructions that are not well-defined for general
    maps between oriented graded posets.
    """

    def __init__(self, source, target, mapping=None,
                 wfcheck=True):
        for x in source, target:
            utils.typecheck(x, {'type': Shape})
        super().__init__(source, target, mapping, wfcheck)

    @classmethod
    def from_ogmap(cls, ogmap):
        return cls(ogmap.source, ogmap.target, ogmap.mapping,
                   wfcheck=False)


def atom(fst, snd):
    return Shape.atom(fst, snd)


def paste(fst, snd, dim):
    return Shape.paste(fst, snd, dim)


def globe(dim):
    return Shape.globe(dim)
