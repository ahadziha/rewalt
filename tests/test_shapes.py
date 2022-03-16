from pytest import raises

from rewal import utils
from rewal.ogposets import (El, OgPoset, OgMap)
from rewal.shapes import (Shape, ShapeMap, atom, paste, globe)


point = Shape.point()
arrow = Shape.arrow()


def test():
    whisker_l = paste(arrow, globe(2))
    whisker_r = paste(globe(2), arrow)
    interch_1 = paste(whisker_l, whisker_r)
    interch_2 = paste(whisker_r, whisker_l)
    interch_3 = paste(globe(2), globe(2), 0)
    assert interch_1 == interch_2 == interch_3


""" Tests for Shape """


def test_Shape():
    assert Shape() == Shape()
    assert Shape().dim == -1


def test_Shape_atom():
    assert Shape.atom(Shape(), Shape()) == point
    assert Shape.atom(point, point).size == [2, 1]
    assert Shape.atom(arrow, arrow).size == [2, 2, 1]


def test_Shape_initial():
    assert Shape.point().initial() == Shape().terminal()
    assert Shape().initial() == Shape().id()


def test_Shape_terminal():
    assert Shape.point().terminal() == Shape.point().id()


""" Tests for ShapeMap """


def test_ShapeMap():
    assert ShapeMap(Shape(), Shape.point()) == \
            ShapeMap(Shape(), Shape.point())
    assert ShapeMap(Shape(), Shape.point()) != \
        OgMap(Shape(), Shape.point())
