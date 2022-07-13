from pytest import raises

from rewalt import utils
from rewalt import shapes
from rewalt.ogposets import (El, OgMap, OgMapPair)
from rewalt.shapes import (Shape, ShapeMap)


""" Tests for Shape """


def test_Shape_init():
    assert Shape() == Shape()
    assert Shape().dim == -1


def test_Shape_isatom():
    empty = Shape.empty()
    arrow = Shape.arrow()
    assert not empty.isatom
    assert arrow.isatom
    assert not arrow.paste(arrow).isatom


def test_Shape_isround():
    frob = Shape.simplex(2).paste(Shape.arrow()).paste(
            Shape.arrow().paste(Shape.simplex(2).dual()))
    assert frob.isround

    whisker_l = Shape.arrow().paste(Shape.globe(2))
    assert not whisker_l.isround


def test_Shape_layers():
    arrow = Shape.arrow()
    globe = Shape.globe(2)
    cospan = globe.paste(arrow).paste(
            arrow.paste(globe), cospan=True)
    shape = cospan.target
    assert shape.layers == [cospan.fst, cospan.snd]


def test_Shape_rewrite_steps():
    arrow = Shape.arrow()
    globe = Shape.globe(2)
    cospan = globe.paste(arrow).paste(
            arrow.paste(globe), cospan=True)
    shape = cospan.target
    assert shape.rewrite_steps == [
            cospan.fst.input,
            cospan.fst.output,
            cospan.snd.output]


def test_Shape_atom():
    empty = Shape.empty()
    point = Shape.point()
    arrow = Shape.arrow()
    globe = Shape.globe(2)
    assert point == empty.atom(empty)
    assert arrow == point.atom(point)

    whisker_l = arrow.paste(globe)
    with raises(ValueError) as err:
        whisker_l.atom(whisker_l)
    assert str(err.value) == utils.value_err(
            whisker_l, 'expecting a round Shape')

    binary = arrow.paste(arrow).atom(arrow)
    cobinary = Shape.simplex(2)
    with raises(ValueError) as err:
        binary.atom(cobinary)
    assert str(err.value) == utils.value_err(
            cobinary, 'input boundary does not match '
            'input boundary of {}'.format(repr(binary)))

    with raises(ValueError) as err:
        cobinary.atom(globe)
    assert str(err.value) == utils.value_err(
            globe, 'output boundary does not match '
            'output boundary of {}'.format(repr(cobinary)))

    assert isinstance(arrow.atom(arrow), shapes.Globe)
    assert isinstance(point.atom(point), shapes.Arrow)
    assert isinstance(binary, shapes.Opetope)

    cospan = arrow.paste(arrow).atom(arrow, cospan=True)
    assert cospan.fst == binary.input
    assert cospan.snd == binary.output


def test_Shape_paste():
    globe = Shape.globe(2)
    arrow = Shape.arrow()

    whisker_l = arrow.paste(globe)
    whisker_r = globe.paste(arrow)

    interch_1 = whisker_l.paste(whisker_r)
    interch_2 = whisker_r.paste(whisker_l)
    interch_3 = globe.paste(globe, 0)
    assert interch_1 == interch_2 == interch_3

    point = Shape.point()
    assert point.paste(arrow, 0) == arrow
    assert globe.paste(arrow, 1) == globe

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

    binary = arrow.paste(arrow).atom(arrow)
    assert isinstance(arrow.paste(arrow), shapes.GlobeString)
    assert isinstance(whisker_l, shapes.Theta)
    assert isinstance(binary.paste(globe), shapes.OpetopeTree)

    cospan = globe.paste(arrow, cospan=True)
    assert cospan.fst.source == globe
    assert cospan.snd.source == arrow
    assert cospan.target == whisker_r


def test_Shape_paste_along():
    arrow = Shape.arrow()
    pasting = arrow.paste(arrow, cospan=True)
    atoming = pasting.target.atom(arrow, cospan=True)
    binary = atoming.target

    fst = binary.output
    snd = pasting.fst.then(atoming.fst)
    snd2 = pasting.snd.then(atoming.fst)

    assoc_l = Shape.paste_along(fst, snd)
    assert isinstance(assoc_l, shapes.OpetopeTree)

    with raises(ValueError) as err:
        Shape.paste_along(snd, snd)
    assert str(err.value) == utils.value_err(
            OgMapPair(snd, snd),
            'not a well-formed span for pasting')

    cospan = Shape.paste_along(fst, snd2, cospan=True)
    assert cospan.fst.image().intersection(
            cospan.snd.image()) == fst.then(cospan.fst).image()


def test_Shape_to_outputs():
    arrow = Shape.arrow()
    twothree = arrow.paste(arrow).atom(arrow.paste(arrow).paste(arrow))
    threetwo = twothree.dual()

    with raises(ValueError) as err:
        twothree.to_outputs([2, 4], twothree)
    assert str(err.value) == utils.value_err(
            [2, 4], 'cannot paste to these outputs')

    with raises(ValueError) as err:
        twothree.to_outputs([2, 3], threetwo)
    assert str(err.value) == utils.value_err(
            [2, 3], 'does not match input boundary of {}'.format(
                repr(threetwo)))

    pasted = twothree.to_outputs([2, 3], twothree)
    assert pasted.size == [7, 8, 2]


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
    assert isinstance(arrow >> point, shapes.Simplex)


def test_Shape_merge():
    arrow = Shape.arrow()
    binary = arrow.paste(arrow).atom(arrow)
    ternary = arrow.paste(arrow).paste(arrow).atom(arrow)
    assoc_l = binary.to_inputs(0, binary)

    assert assoc_l.merge() == ternary
    assert isinstance(assoc_l.merge(), shapes.Opetope)


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
