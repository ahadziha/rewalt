from pytest import raises

from rewal import utils
from rewal.ogposets import (OgMap)
from rewal.shapes import (Shape, ShapeMap, atom, paste, globe)


empty = Shape()
point = globe(0)
arrow = globe(1)
globe2 = globe(2)

whisker_l = paste(arrow, globe2)
whisker_r = paste(globe2, arrow)

binary = atom(paste(arrow, arrow), arrow)
cobinary = atom(arrow, paste(arrow, arrow))

assoc_l = paste(paste(binary, arrow), binary)
assoc_r = paste(paste(arrow, binary), binary)
associator = atom(assoc_l, assoc_r)

frob_l = paste(paste(cobinary, arrow), paste(arrow, binary))
frob_r = paste(paste(arrow, cobinary), paste(binary, arrow))
frob_c = paste(binary, cobinary)
frobenius1 = atom(frob_l, frob_r)
frobenius2 = atom(frob_c, frob_l)
frobenius3 = atom(frob_r, frob_c)


def test():
    assert associator.size == frobenius1.size == [4, 6, 4, 1]


""" Tests for Shape """


def test_Shape():
    assert Shape() == Shape()
    assert Shape().dim == -1


def test_Shape_isround():
    assert cobinary.isround
    assert frob_l.isround
    assert not whisker_l.isround


def test_Shape_atom():
    assert atom(empty, empty) == point
    assert atom(point, point) == arrow
    assert atom(arrow, arrow) == globe2

    with raises(ValueError) as err:
        atom(whisker_l, whisker_l)
    assert str(err.value) == utils.value_err(
            whisker_l, 'expecting a round Shape')

    with raises(ValueError) as err:
        atom(binary, cobinary)
    assert str(err.value) == utils.value_err(
            cobinary, 'input boundary does not match '
            'input boundary of {}'.format(repr(binary)))

    with raises(ValueError) as err:
        atom(cobinary, globe2)
    assert str(err.value) == utils.value_err(
            globe2, 'output boundary does not match '
            'output boundary of {}'.format(repr(cobinary)))


def test_Shape_paste():
    interch_1 = paste(whisker_l, whisker_r)
    interch_2 = paste(whisker_r, whisker_l)
    interch_3 = paste(globe2, globe2, 0)
    assert interch_1 == interch_2 == interch_3

    assert paste(arrow, point, 0) == arrow == paste(point, arrow, 0)
    assert paste(globe2, arrow, 1) == globe2 == paste(arrow, globe2, 1)

    with raises(ValueError) as err:
        paste(point, point)
    assert str(err.value) == utils.value_err(
            -1, 'expecting non-negative integer')

    with raises(ValueError) as err:
        paste(binary, binary)
    assert str(err.value) == utils.value_err(
            binary,
            'input 1-boundary does not match '
            'output 1-boundary of {}'.format(repr(binary)))


def test_Shape_suspend():
    assert Shape.suspend(whisker_l).size == [2] + whisker_l.size
    assert Shape.suspend(arrow) == globe2


def test_Shape_initial():
    assert point.initial() == empty.terminal()
    assert empty.initial() == empty.id()
    assert whisker_l.initial().istotal


def test_Shape_terminal():
    assert point.terminal() == point.id()
    assert whisker_l.terminal().issurjective


""" Tests for ShapeMap """


def test_ShapeMap_init():
    undefined = OgMap(point, arrow)
    with raises(ValueError) as err:
        ShapeMap(undefined)
    assert str(err.value) == utils.value_err(
            undefined,
            'a ShapeMap must be total')
