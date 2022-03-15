"""
Implements shapes of cells and diagrams.
"""

from rewal import utils
from rewal.ogposets import (El, OgPoset, GrSet, GrSubset, OgMap)


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
        boundary_in, boundary_out = self.boundary('-'), self.boundary('+')
        intersection = boundary_in.intersection(boundary_out)
        for k in range(self.dim-2, -1, -1):
            boundary_in = boundary_in.boundary('-')
            boundary_out = boundary_out.boundary('+')
            if not intersection.issubset(boundary_in.union(boundary_out)):
                return False
            intersection = boundary_in.intersection(boundary_out)
        return True

    @classmethod
    def point(cls):
        """
        Generates the point shape.
        """
        point = cls.__new__(cls)
        point._face_data = [
                [{'-': set(), '+': set()}]
                ]
        point._coface_data = [
                [{'-': set(), '+': set()}]
                ]
        point._isatom = True
        return point
    
    def boundary_inclusion(self, sign=None, dim=None):
        """
        Input and output boundaries of Shapes are Shapes.
        """
        boundary_ogmap = super().boundary_inclusion(sign, dim)
        if sign is None:
            return boundary_ogmap
        reordering = self.__reorder(boundary_ogmap.source)
        return ShapeMap.from_ogmap(
                reordering.then(boundary_ogmap))

    def id(self):
        return ShapeMap.from_ogmap(super().id())

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
                                x for x in mapping[focus.dim - 1]
                                if len(candidates(x)) > 0)
                        focus_stack.append(GrSubset(
                            GrSet(candidates(x)[0]),
                            shape, wfcheck=False).closure())

        def reordered_faces(x, sign):
            return {k for k in range(shape.size[x.dim - 1])
                    if mapping[x.dim - 1][k] in shape.faces(x, sign)}

        def reordered_cofaces(x, sign):
            return {k for k in range(shape.size[x.dim + 1])
                    if mapping[x.dim + 1][k] in shape.cofaces(x, sign)}

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
