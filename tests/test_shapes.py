import numpy as np
from pytest import raises

from rewal import utils
from rewal.shapes import El, OgPos


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


# Test on the "right-whiskered 2-globe".

test_face = [
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
test_coface = [
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

test_ogpos = OgPos(test_face, test_coface)


def test_OgPos_size():
    assert test_ogpos.size == [3, 3, 1]


def test_OgPos_dim():
    assert test_ogpos.dim == 2


def test_OgPos_chain():
    chain = [
            np.array([[-1, -1, 0], [1, 1, -1], [0, 0, 1]]),
            np.array([[-1], [1], [0]])
            ]
    test_chain = test_ogpos.chain
    assert (test_chain[0] == chain[0]).all() and \
        (test_chain[1] == chain[1]).all()


def test_OgPos_getitem():
    assert test_ogpos[1] == [0, 1, 2]

    with raises(TypeError) as err:
        test_ogpos['x']
    assert str(err.value) == utils.type_err(int, 'x')

    with raises(ValueError) as err:
        test_ogpos[3]
    assert str(err.value) == utils.value_err(3, 'out of bounds')


def test_OgPos_from_face_data():
    assert test_ogpos == OgPos.from_face_data(test_face)
