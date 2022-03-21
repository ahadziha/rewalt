from pytest import raises

from rewal import utils
from rewal.ogposets import (El, OgMap)
from rewal.shapes import (Shape, ShapeMap)


""" Tests for Shape """


def test_Shape_init():
    assert Shape() == Shape()
    assert Shape().dim == -1


def test_Shape_isround():
    frob = Shape.simplex(2).paste(Shape.arrow()).paste(
            Shape.arrow().paste(Shape.simplex(2).dual()))
    assert frob.isround

    whisker_l = Shape.arrow().paste(Shape.globe(2))
    assert not whisker_l.isround


def test_Shape_atom():
    point = Shape.point()
    assert point == Shape.empty().atom(Shape.empty())
    assert Shape.arrow() == point.atom(point)

    whisker_l = Shape.arrow().paste(Shape.globe(2))
    with raises(ValueError) as err:
        whisker_l.atom(whisker_l)
    assert str(err.value) == utils.value_err(
            whisker_l, 'expecting a round Shape')

    binary = Shape.simplex(2).dual()
    cobinary = Shape.simplex(2)
    with raises(ValueError) as err:
        binary.atom(cobinary)
    assert str(err.value) == utils.value_err(
            cobinary, 'input boundary does not match '
            'input boundary of {}'.format(repr(binary)))

    globe2 = Shape.globe(2)
    with raises(ValueError) as err:
        cobinary.atom(globe2)
    assert str(err.value) == utils.value_err(
            globe2, 'output boundary does not match '
            'output boundary of {}'.format(repr(cobinary)))


def test_Shape_paste():
    globe2 = Shape.globe(2)
    arrow = Shape.arrow()

    whisker_l = arrow.paste(globe2)
    whisker_r = globe2.paste(arrow)

    interch_1 = whisker_l.paste(whisker_r)
    interch_2 = whisker_r.paste(whisker_l)
    interch_3 = globe2.paste(globe2, 0)
    assert interch_1 == interch_2 == interch_3

    point = Shape.point()
    assert point.paste(arrow, 0) == arrow
    assert globe2.paste(arrow, 1) == globe2

    with raises(ValueError) as err:
        point.paste(point)
    assert str(err.value) == utils.value_err(
            -1, 'expecting non-negative integer')

    cobinary = Shape.simplex(2)
    with raises(ValueError) as err:
        cobinary.paste(cobinary)
    assert str(err.value) == utils.value_err(
            cobinary,
            'input 1-boundary does not match '
            'output 1-boundary of {}'.format(repr(cobinary)))


def test_Shape_suspend():
    globe2 = Shape.globe(2)
    arrow = Shape.arrow()

    whisker_l = arrow.paste(globe2)

    assert whisker_l.suspend().size == [2] + whisker_l.size
    assert arrow.suspend(2) == globe2.suspend()


def test_Shape_gray():
    globe2 = Shape.globe(2)
    arrow = Shape.arrow()

    assert (globe2 * arrow).size == [4, 6, 4, 1]
    assert Shape.gray() == Shape.point()
    assert (arrow * arrow) * arrow == arrow * (arrow * arrow)


def test_Shape_inflate():
    simplex2 = Shape.simplex(2)
    assert simplex2.inflate().source.boundary_inclusion('-').source == \
        simplex2


def test_Shape_under():
    globe2 = Shape.globe(2)
    arrow = Shape.arrow()
    whisker_l = arrow.paste(globe2)

    assert whisker_l.under(El(2, 0)).source == globe2


def test_Shape_initial():
    point = Shape.point()
    empty = Shape.empty()

    assert point.initial() == empty.terminal()
    assert empty.initial() == empty.id()
    assert point.initial().istotal


def test_Shape_terminal():
    point = Shape.point()
    assert point.terminal() == point.id()


""" Tests for ShapeMap """


def test_ShapeMap_init():
    point = Shape.point()
    arrow = Shape.arrow()

    undefined = OgMap(point, arrow)
    with raises(ValueError) as err:
        ShapeMap(undefined)
    assert str(err.value) == utils.value_err(
            undefined,
            'a ShapeMap must be total')
