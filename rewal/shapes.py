"""
Implements shapes of cells and diagrams.
"""

from rewal import utils
from rewal.ogposets import OgPoset


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

    @classmethod
    def __upgrade(cls, ogposet):
        """
        Private method 'forcing' upgrade of an OgPoset to a shape class.
        """
        utils.typecheck(ogposet, {'type': OgPoset})
        shape = cls.__new__(cls)
        shape._face_data = ogposet.face_data
        shape._coface_data = ogposet.coface_data
        return shape
