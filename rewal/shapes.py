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
        super().__init__([], [], wfcheck=False)

    # Redefining to be more lax wrt subclasses
    def __eq__(self, other):
        return isinstance(other, Shape) and \
                self.face_data == other.face_data and \
                self.coface_data == other.coface_data

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
                OgPoset(face_data, coface_data, wfcheck=False)).source
        return new_atom

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
            return OgMapPair(span.snd, snd.id())
        if dim >= snd.dim:
            return OgMapPair(fst.id(), span.fst)
        pushout = span.pushout()
        reordering = Shape.__reorder(pushout.target).inv()
        paste_cospan = pushout.then(reordering)

        # Theta is closed under pasting
        if isinstance(fst, Theta) and isinstance(snd, Theta):
            target = Theta.__upgrade(paste_cospan.target)
            paste_cospan = OgMapPair(
                    OgMap(fst, target, paste_cospan.fst.mapping,
                          wfcheck=False),
                    OgMap(snd, target, paste_cospan.snd.mapping,
                          wfcheck=False))

        return OgMapPair(
                ShapeMap(paste_cospan.fst, wfcheck=False),
                ShapeMap(paste_cospan.snd, wfcheck=False))

    @staticmethod
    def paste(fst, snd, dim=None):
        return Shape.paste_cospan(fst, snd, dim).target

    # Other constructors
    def suspend(self, n=1):
        if n == 0:
            return self

        return Shape.__reorder(
                OgPoset.suspension(self, n)).source

    # Special maps
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
        return ShapeMap(reordering.then(boundary_ogmap),
                        wfcheck=False)

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
    @classmethod
    def point(cls):
        """ The point. """
        face_data = [[{'-': set(), '+': set()}]]
        coface_data = [[{'-': set(), '+': set()}]]

        return cls.__upgrade(
                OgPoset(face_data, coface_data,
                        wfcheck=False, matchcheck=False))

    @classmethod
    def arrow(cls):
        """ The arrow. """
        return cls.__upgrade(Globe(1))

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
                dim = focus.support.dim
                if dim == 0:
                    for x in focus[dim]:
                        mapping[dim].append(x)
                        marked.support.add(x)
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
                                        if y not in marked]
                            x = next(
                                    x for x in mapping[dim - 1]
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
    def __upgrade(cls, ogp):
        """
        Forces upgrade of an OgPoset to the shape class.
        """
        shape = cls.__new__(cls)
        shape._face_data = ogp.face_data
        shape._coface_data = ogp.coface_data
        return shape


class Theta(Shape):
    """
    Class for Batanin cell shapes (objects of the Theta category).
    """
    def __init__(self, *thetas):
        def tree(*thetas):
            if thetas:
                theta = thetas[0]
                utils.typecheck(theta, {'type': Theta})
                thetas = thetas[1:]
                if thetas:
                    return Shape.paste(
                            theta.suspend(),
                            tree(*thetas), 0)
                return theta.suspend()
            return Shape.point()
        new = tree(*thetas)

        super().__init__()
        self._face_data = new.face_data
        self._coface_data = new.coface_data

    # Theta and Globe are closed under suspension
    def suspend(self, dim=1):
        return self.__class__._Shape__upgrade(
                super().suspend(dim))


class Globe(Theta):
    """ The globes. """
    def __init__(self, dim):
        utils.typecheck(dim, {
            'type': int,
            'st': lambda n: n >= 0,
            'why': 'expecting non-negative integer'})
        new = Shape.point().suspend(dim)

        super().__init__()
        self._face_data = new.face_data
        self._coface_data = new.coface_data


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


def atom(fst, snd):
    return Shape.atom(fst, snd)


def paste(fst, snd, dim=None):
    return Shape.paste(fst, snd, dim)
