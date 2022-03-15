"""
Implements shapes of cells and diagrams.
"""

from rewal import utils
from rewal.ogposets import (OgPoset, OgMap)


class Shape(OgPoset):
    """
    Class for shapes of cells and diagrams.
    """

    def __init__(self):
        self._face_data = []
        self._coface_data = []

    # Class methods
    @classmethod
    def point(cls):
        """
        Generates the point shape.
        This is a class method because the point also belongs to all
        specialised shape classes.
        """
        point = cls.__new__(cls)
        point._face_data = [
                [{'-': set(), '+': set()}]
                ]
        point._coface_data = [
                [{'-': set(), '+': set()}]
                ]
        return point

    # Internal methods
    @classmethod
    def _upgrade(cls, ogposet):
        """
        Forces upgrade of an OgPoset to the shape class.
        """
        utils.typecheck(ogposet, {'type': OgPoset})
        shape = cls.__new__(cls)
        shape._face_data = ogposet.face_data
        shape._coface_data = ogposet.coface_data
        return shape


class ShapeMap(OgMap):
    """
    An OgMap overlay used for maps between Shape objects.
    Used for constructions that are not well-defined for general
    maps between oriented graded posets.
    """

    def __init__(self, source, target, mapping=None,
                 wfcheck=True):
        for x in source, target:
            utils.typecheck(x, {'type': Shape})
        super().__init__(source, target, mapping, wfcheck)

    # Internal methods
    @staticmethod
    def _upgrade(ogmap):
        """
        Forces upgrade of an OgMap to a ShapeMap.
        """
        utils.typecheck(ogmap, {'type': OgMap})
        return ShapeMap(
                Shape._upgrade(ogmap.source),
                Shape._upgrade(ogmap.target),
                ogmap.mapping,
                wfcheck=False)
