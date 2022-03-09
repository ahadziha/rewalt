import numpy as np
from pytest import raises

from rewal import utils
from rewal.shapes import El, OgPos, OgMap


# El tests

def test_El_init():
    assert El(2, 3) == El(2, 3)

    with raises(TypeError) as err:
        El('x', 2)
    assert str(err.value) == utils.type_err(int, 'x')

    with raises(ValueError) as err:
        El(3, -1)
    assert str(err.value) == utils.value_err(
            -1, 'expecting non-negative integer')


def test_El_shift():
    assert El(2, 3).shift(4) == El(2, 7)

    with raises(TypeError) as err:
        El(2, 3).shift('x')
    assert str(err.value) == utils.type_err(int, 'x')

    with raises(ValueError) as err:
        El(2, 3).shift(-4)
    assert str(err.value) == utils.value_err(-4, 'out of bounds')


# OgPos tests

def test_OgPos_init():
    test_face = [
        [
            {'-': set(), '+': set()},
            {'-': set(), '+': set()},
        ], [
            {'-': {0}, '+': {1}}
        ]]
    test_coface = [
        [
            {'-': {0}, '+': set()},
            {'-': set(), '+': {0}},
        ], [
            {'-': set(), '+': set()}
        ]]

    assert OgPos(test_face, test_coface) == OgPos(test_face, test_coface)

    test_face[0][0]['-'] = {0}
    with raises(ValueError) as err:
        OgPos(test_face, test_coface)
    assert str(err.value) == utils.value_err(0, 'out of bounds')

    test_face[0][0]['-'] = set()
    test_face[1][0]['-'] = {2}
    with raises(ValueError) as err:
        OgPos(test_face, test_coface)
    assert str(err.value) == utils.value_err(2, 'out of bounds')

    test_face[1][0]['-'] = {1}
    with raises(ValueError) as err:
        OgPos(test_face, test_coface)
    assert str(err.value) == 'Input and output faces must be disjoint.'

    test_face[1][0]['-'] = set()
    test_face[1][0]['+'] = set()
    with raises(ValueError) as err:
        OgPos(test_face, test_coface)
    assert str(err.value) == 'The element must have at least one face.'

    test_face[1][0]['-'] = {1}
    test_face[1][0]['+'] = {0}
    with raises(ValueError) as err:
        OgPos(test_face, test_coface)
    assert str(err.value) == 'Face and coface data do not match.'

    test_face[1][0]['-'] = {'x'}
    with raises(TypeError) as err:
        OgPos(test_face, test_coface)
    assert str(err.value) == utils.type_err(int, 'x')

    test_face[1][0]['-'] = 0
    with raises(TypeError) as err:
        OgPos(test_face, test_coface)
    assert str(err.value) == utils.type_err(set, 0)

    test_face[1][0] = {'k': {0}}
    with raises(ValueError) as err:
        OgPos(test_face, test_coface)
    assert str(err.value) == utils.value_err(
            {'k': {0}},
            "expecting dict with keys '-', '+'")

    test_face[1][0] = {0}
    with raises(TypeError) as err:
        OgPos(test_face, test_coface)
    assert str(err.value) == utils.type_err(dict, {0})

    test_face[1] = []
    with raises(ValueError) as err:
        OgPos(test_face, test_coface)
    assert str(err.value) == utils.value_err(
            [], 'expecting non-empty list')

    test_face[1] = {0}
    with raises(TypeError) as err:
        OgPos(test_face, test_coface)
    assert str(err.value) == utils.type_err(list, {0})

    test_face = {0}
    with raises(TypeError) as err:
        OgPos(test_face, test_coface)
    assert str(err.value) == utils.type_err(list, {0})


# The next tests will be done on the "right-whiskered 2-globe".

whisker_face = [
        [
            {'-': set(), '+': set()},
            {'-': set(), '+': set()},
            {'-': set(), '+': set()}
        ], [
            {'-': {0}, '+': {1}},
            {'-': {0}, '+': {1}},
            {'-': {1}, '+': {2}}
        ], [
            {'-': {0}, '+': {1}}
        ]]
whisker_coface = [
        [
            {'-': {0, 1}, '+': set()},
            {'-': {2}, '+': {0, 1}},
            {'-': set(), '+': {2}}
        ], [
            {'-': {0}, '+': set()},
            {'-': set(), '+': {0}},
            {'-': set(), '+': set()}
        ], [
            {'-': set(), '+': set()}
        ]]

whisker = OgPos(whisker_face, whisker_coface)


def test_OgPos_size():
    assert whisker.size == [3, 3, 1]


def test_OgPos_dim():
    assert whisker.dim == 2


def test_OgPos_chain():
    chain = [
            np.array([[-1, -1, 0], [1, 1, -1], [0, 0, 1]]),
            np.array([[-1], [1], [0]])
            ]
    test_chain = whisker.chain
    assert (test_chain[0] == chain[0]).all() and \
        (test_chain[1] == chain[1]).all()


def test_OgPos_getitem():
    assert whisker[1] == [0, 1, 2]

    with raises(TypeError) as err:
        whisker['x']
    assert str(err.value) == utils.type_err(int, 'x')

    with raises(ValueError) as err:
        whisker[3]
    assert str(err.value) == utils.value_err(3, 'out of bounds')


def test_OgPos_from_face_data():
    assert whisker == OgPos.from_face_data(whisker_face)


# Introducing the interval (1-globe).

interval_face = [
        [
            {'-': set(), '+': set()},
            {'-': set(), '+': set()},
        ], [
            {'-': {0}, '+': {1}},
        ]]

interval = OgPos.from_face_data(interval_face)


# Tests on OgMap

def test_OgMap_init():
    assert OgMap(whisker, interval) == OgMap(whisker, interval)


def test_OgMap_mapping():
    assert OgMap(interval, whisker).mapping == [[None, None], [None]]
