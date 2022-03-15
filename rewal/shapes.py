"""
Implements shapes of cells and diagrams.
"""

from rewal import utils
from rewal.ogposets import (El, OgPoset, OgMap)


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
        """
        return self._isatom

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
    @classmethod
    def __upgrade(cls, ogposet):
        """
        Forces upgrade of an OgPoset to the shape class.
        """
        shape = cls.__new__(cls)
        shape._face_data = ogposet.face_data
        shape._coface_data = ogposet.coface_data
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

    @staticmethod
    def from_ogmap(ogmap):
        return ShapeMap(ogmap.source, ogmap.target, ogmap.mapping,
                        wfcheck=False)
