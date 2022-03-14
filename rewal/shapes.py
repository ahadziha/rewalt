"""
Implements shapes of cells and diagrams.
"""

from rewal.ogposets import OgPoset


class Shape(OgPoset):
    """
    Class for shapes of cells and diagrams.
    """

    def __init__(self):
        self._face_data = [
                [{'-': set(), '+': set()}]
                ]
        self._coface_data = [
                [{'-': set(), '+': set()}]
                ]

    @staticmethod
    def point():
        return Shape()

    # Private methods
    @staticmethod
    def __force(ogposet):
        """
        Private method to 'force' upgrade of an OgPoset to a Shape.
        """
        shape = Shape()
        shape._face_data = ogposet.face_data
        shape._coface_data = ogposet.coface_data

        return shape
