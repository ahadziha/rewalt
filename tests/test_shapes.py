from pytest import raises

from rewal import utils
from rewal.ogposets import (El, OgMap)
from rewal.shapes import (Shape, ShapeMap, Simplex, Opetope)


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


def test_Shape_join():
    point = Shape.point()
    arrow = Shape.arrow()

    assert arrow >> point == point << point << point
    assert Shape.join() == Shape.empty()
    assert isinstance(arrow >> point, Simplex)


def test_Shape_merge():
    arrow = Shape.arrow()
    binary = arrow.paste(arrow).atom(arrow)
    ternary = arrow.paste(arrow).paste(arrow).atom(arrow)
    assoc_l = binary.to_inputs(0, binary)

    assert assoc_l.merge() == ternary
    assert isinstance(assoc_l.merge(), Opetope)


def test_Shape_inflate():
    simplex2 = Shape.simplex(2)
    assert simplex2.inflate().source.boundary('-').source == \
        simplex2

    whisker_l = Shape.arrow().paste(Shape.globe(2))
    assert whisker_l.inflate().source.boundary('+').source == \
        whisker_l


def test_Shape_atom_inclusion():
    globe2 = Shape.globe(2)
    arrow = Shape.arrow()
    whisker_l = arrow.paste(globe2)

    assert whisker_l.atom_inclusion(El(2, 0)).source == globe2


def test_Shape_initial():
    point = Shape.point()
    empty = Shape.empty()

    assert point.initial() == empty.terminal()
    assert empty.initial() == empty.id()
    assert point.initial().istotal


def test_Shape_terminal():
    point = Shape.point()
    assert point.terminal() == point.id()


""" Tests for Shape subclasses """


def test_Simplex():
    arrow = Shape.simplex(1)
    triangle = Shape.simplex(2)
    tetra = Shape.simplex(3)

    assert tetra.simplex_face(3).source == triangle
    assert tetra.simplex_face(0) == tetra.atom_inclusion(
            El(2, 3))

    map1 = triangle.simplex_degeneracy(2).then(
        arrow.simplex_degeneracy(1))
    map2 = triangle.simplex_degeneracy(1).then(
        arrow.simplex_degeneracy(1))
    assert map1 == map2

    map3 = tetra.simplex_face(2).then(
        triangle.simplex_degeneracy(2))
    assert map3 == triangle.id()

    map4 = tetra.simplex_face(0).then(
        triangle.simplex_degeneracy(2))
    map5 = arrow.simplex_degeneracy(1).then(
        triangle.simplex_face(0))
    assert map4 == map5


def test_Cube():
    arrow = Shape.cube(1)
    square = Shape.cube(2)
    cube = Shape.cube(3)

    map1 = square.cube_degeneracy(2).then(
        arrow.cube_degeneracy(1))
    map2 = square.cube_degeneracy(1).then(
        arrow.cube_degeneracy(1))
    assert map1 == map2

    map3 = square.cube_face(0, '+').then(
        cube.cube_face(2, '-'))
    map4 = square.cube_face(1, '-').then(
        cube.cube_face(0, '+'))
    assert map3 == map4

    map5 = square.cube_connection(1, '-').then(
        arrow.cube_connection(0, '-'))
    map6 = square.cube_connection(0, '-').then(
        arrow.cube_connection(0, '-'))
    assert map5 == map6


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
