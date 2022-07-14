from pytest import raises

from rewalt import utils, diagrams
from rewalt.shapes import Shape
from rewalt.diagrams import DiagSet


""" Tests for DiagSet """


def test_DiagSet_init():
    X = DiagSet()
    assert X == X
    Y = DiagSet()
    assert not X == Y
    assert len(X) == 0


""" Some example objects """

Mon = DiagSet()
pt = Mon.add('pt')
a = Mon.add('a', pt, pt)
m = Mon.add('m', a.paste(a), a)
u = Mon.add('u', pt.unit(), a)
assoc = Mon.add(
        'assoc',
        m.to_inputs(0, m),
        m.to_inputs(1, m))
lunit = Mon.add(
        'lunit',
        m.to_inputs(0, u),
        a.lunitor('-'))
runit = Mon.add(
        'runit',
        m.to_inputs(1, u),
        a.runitor('-'))

RP3 = DiagSet()
c0 = RP3.add_simplex('c0')
c1 = RP3.add_simplex('c1', c0, c0)
c2 = RP3.add_simplex(
        'c2',
        c1,
        c0.simplex_degeneracy(0),
        c1)
c3 = RP3.add_simplex(
        'c3',
        c2,
        c1.simplex_degeneracy(0),
        c1.simplex_degeneracy(1),
        c2)

T2 = DiagSet()
v = T2.add_cube('v')
e0 = T2.add_cube('e0', v, v)
e1 = T2.add_cube('e1', v, v)
s = T2.add_cube('s', e0, e0, e1, e1)

C = DiagSet()
x = C.add('x')
y = C.add('y')
z = C.add('z')
f = C.add('f', x, y)
g = C.add('g', y, z)
h, c_fg = C.compose(f.paste(g), 'h', 'c_fg')


def test_DiagSet_str():
    assert str(Mon) == 'DiagSet with 7 generators'


def test_DiagSet_getitem():
    assert m == Mon['m']
    assert assoc == Mon['assoc']
    with raises(KeyError) as err:
        Mon['n']
    assert str(err.value) == "'n'"


def test_DiagSet_contains():
    assert 'm' in Mon
    assert 'n' not in Mon


def test_DiagSet_len():
    assert len(Mon) == 7


def test_DiagSet_generators():
    X = DiagSet()
    X.add('x')
    assert X.generators == {
            'x': {
                'shape': Shape.point(),
                'mapping': [['x']],
                'faces': set(),
                'cofaces': set()}}


def test_DiagSet_by_dim():
    assert Mon.by_dim == {
            0: {'pt'},
            1: {'a'},
            2: {'m', 'u'},
            3: {'assoc', 'lunit', 'runit'}}


def test_DiagSet_compositors():
    assert C.compositors == {
            'c_fg': {
                'shape': f.paste(g).shape,
                'mapping': f.paste(g).mapping}}


def test_DiagSet_dim():
    assert Mon.dim == 3
    assert C.dim == 2


def test_DiagSet_issimplicial():
    assert RP3.issimplicial
    assert not C.issimplicial


def test_DiagSet_iscubical():
    assert T2.iscubical
    assert not Mon.iscubical


def test_DiagSet_add():
    X = DiagSet()
    x = X.add('x')

    with raises(ValueError) as err:
        X.add('x')
    assert str(err.value) == utils.value_err(
            'x', 'name already in use')

    y = X.add('y', linvertor='test')
    with raises(KeyError) as err:
        X.generators['y']['linvertor']
    assert str(err.value) == "'linvertor'"

    with raises(ValueError) as err:
        X.add('a', pt, pt)
    assert str(err.value) == utils.value_err(
            pt, 'not a diagram in {}'.format(repr(X)))

    a = X.add('a', x, y)
    b = X.add('b', y, x)
    with raises(ValueError) as err:
        X.add('p', a, b)
    assert str(err.value) == utils.value_err(
            b,
            'boundary does not match boundary of {}'.format(repr(a)))


def test_DiagSet_add_simplex():
    S = DiagSet()
    x0 = S.add_simplex('x0')

    with raises(ValueError) as err:
        S.add_simplex('x0')
    assert str(err.value) == utils.value_err(
            'x0', 'name already in use')

    x1 = S.add_simplex('x1', x0, x0)
    paste = x1.paste(x1)
    with raises(TypeError) as err:
        S.add_simplex('x2', paste, x1)
    assert str(err.value) == utils.type_err(
            diagrams.SimplexDiagram, paste)

    with raises(ValueError) as err:
        S.add_simplex('x2', x1, x0, x1)
    assert str(err.value) == utils.value_err(
            x0, 'expecting a 1-simplex in {}'.format(
                repr(S)))

    y0 = S.add_simplex('y0')
    y1 = S.add_simplex('y1', y0, x0)
    with raises(ValueError) as err:
        S.add_simplex('x2', x1, y1, y1)
    assert str(err.value) == utils.value_err(
            y1, 'boundary of face does not match other faces')
