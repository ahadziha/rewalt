from pytest import raises

from rewalt import utils, diagrams
from rewalt.ogposets import El
from rewalt.shapes import Shape
from rewalt.diagrams import (DiagSet, Diagram)


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
assinv, rinv, linv = Mon.invert('assoc')

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
    assert str(Mon) == 'DiagSet with 10 generators'


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
    assert len(Mon) == 10


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
            3: {'assoc', 'lunit', 'runit', 'assoc⁻¹'},
            4: {'inv(assoc, assoc⁻¹)', 'inv(assoc⁻¹, assoc)'}}


def test_DiagSet_compositors():
    assert C.compositors == {
            'c_fg': {
                'shape': f.paste(g).shape,
                'mapping': f.paste(g).mapping}}


def test_DiagSet_dim():
    assert Mon.dim == 4
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


def test_DiagSet_add_cube():
    K = DiagSet()
    x0 = K.add_cube('x0')

    with raises(ValueError) as err:
        K.add_cube('x0')
    assert str(err.value) == utils.value_err(
            'x0', 'name already in use')

    x1 = K.add_cube('x1', x0, x0)
    paste = x1.paste(x1)
    with raises(TypeError) as err:
        K.add_cube('x2', paste, paste)
    assert str(err.value) == utils.type_err(
            diagrams.CubeDiagram, paste)

    with raises(ValueError) as err:
        K.add_cube('x2', x1, x1, x0, x1)
    assert str(err.value) == utils.value_err(
            x0, 'expecting a 1-cube in {}'.format(
                repr(K)))

    y0 = K.add_cube('y0')
    y1 = K.add_cube('y1', x0, y0)
    with raises(ValueError) as err:
        K.add_cube('x2', x1, x1, y1, y1)
    assert str(err.value) == utils.value_err(
            y1, 'boundary of face does not match other faces')


def test_DiagSet_invert():
    X = DiagSet()
    x = X.add('x')
    y = X.add('y')
    f = X.add('f', x, y)
    g, rinv, linv = X.invert('f')

    assert g.input == y
    assert g.output == x
    assert rinv.input == f.paste(g)
    assert rinv.output == x.unit()
    assert linv.input == g.paste(f)
    assert linv.output == y.unit()
    assert g.name == 'f⁻¹'
    assert rinv.name == 'inv(f, f⁻¹)'
    assert linv.name == 'inv(f⁻¹, f)'

    with raises(ValueError) as err:
        X.invert('x')
    assert str(err.value) == utils.value_err(
            'x', 'cannot invert 0-cell')

    with raises(ValueError) as err:
        X.invert('f⁻¹')
    assert str(err.value) == utils.value_err(
            'f⁻¹', 'already inverted')


def test_DiagSet_make_inverses():
    X = DiagSet()
    x = X.add('x')
    y = X.add('y')
    f = X.add('f', x, y)
    g = X.add('g', y, x)
    rinv, linv = X.make_inverses('f', 'g')

    assert rinv.input == f.paste(g)
    assert rinv.output == x.unit()
    assert linv.input == g.paste(f)
    assert linv.output == y.unit()
    assert rinv.name == 'inv(f, g)'
    assert linv.name == 'inv(g, f)'

    X.add('h', y, x)
    with raises(ValueError) as err:
        X.make_inverses('h', 'f')
    assert str(err.value) == utils.value_err(
            'f', 'already inverted')


def test_DiagSet_compose():
    X = DiagSet()
    x = X.add('x')
    y = X.add('y')
    z = X.add('z')
    f = X.add('f', x, y)
    g = X.add('g', y, z)

    fpasteg = f.paste(g)
    fg, c_fg = X.compose(f.paste(g))
    assert fg.input == x
    assert fg.output == z
    assert c_fg.input == f.paste(g)
    assert c_fg.output == fg
    assert fg.name == '⟨{}⟩'.format(f.paste(g).name)
    assert c_fg.name == 'comp({})'.format(f.paste(g).name)

    p = X.add('p', f, f)
    q = X.add('q', g, g)

    notround = p.paste(q, 0)
    with raises(ValueError) as err:
        X.compose(notround)
    assert str(err.value) == utils.value_err(
            notround, 'composable diagrams must have round shape')

    with raises(ValueError) as err:
        X.compose(fpasteg)
    assert str(err.value) == utils.value_err(
            fpasteg, 'already has a composite')


def test_DiagSet_make_composite():
    X = DiagSet()
    x = X.add('x')
    y = X.add('y')
    z = X.add('z')
    f = X.add('f', x, y)
    g = X.add('g', y, z)
    h = X.add('h', x, z)

    c_fg = X.make_composite('h', f.paste(g))
    assert c_fg.input == f.paste(g)
    assert c_fg.output == h
    assert c_fg.name == 'comp({})'.format(f.paste(g).name)

    p = X.add('p', f, f)
    q = X.add('q', g, g)
    X.add('r', f.paste(g), f.paste(g))

    notround = p.paste(q, 0)
    with raises(ValueError) as err:
        X.make_composite('r', notround)
    assert str(err.value) == utils.value_err(
            notround, 'composable diagrams must have round shape')

    X.add('k', x, z)
    fpasteg = f.paste(g)
    with raises(ValueError) as err:
        X.make_composite('k', fpasteg)
    assert str(err.value) == utils.value_err(
            fpasteg, 'already has a composite')


def test_DiagSet_remove():
    X = DiagSet()
    x = X.add('x')
    y = X.add('y')
    X.add('f', x, y)
    assert 'y' in X
    assert 'f' in X

    X.remove('y')
    assert 'y' not in X
    assert 'f' not in X
    assert 'x' in X

    a = X.add('a', x, x)
    X.add('b', x, x)
    rinv, linv = X.make_inverses('a', 'b')
    assert a.isinvertiblecell
    X.remove('b')
    assert not a.isinvertiblecell


def test_DiagSet_update():
    X = DiagSet()
    X.add('x', color='blue')
    assert X.generators['x']['color'] == 'blue'
    X.update('x', color='magenta')
    assert X.generators['x']['color'] == 'magenta'

    with raises(AttributeError) as err:
        X.update('x', shape='circle')
    assert str(err.value) == "('shape', 'private attribute')"


def test_DiagSet_yoneda():
    arrow = Shape.arrow()
    embedarrow = DiagSet.yoneda(arrow)

    assert embedarrow.by_dim == {
            0: {El(0, 0), El(0, 1)},
            1: {El(1, 0)}}


""" Tests for Diagram """


def test_Diagram_init():
    empty = Diagram(C)
    assert empty.shape == Shape.empty()
    assert empty.mapping == []
    assert empty.ambient == C


def test_Diagram_str():
    assert str(a) == 'a'
    assert str(c1) == 'c1'


def test_Diagram_eq():
    assert a == Mon['a']
    assert c2.output == c1.paste(c1)


def test_Diagram_len():
    assert len(a) == 3


def test_Diagram_getitem():
    assert c2[El(2, 0)] == 'c2'
    assert c2[El(1, 0)] == 'c0'
    assert c2[El(1, 1)] == 'c1'


def test_Diagram_contains():
    assert El(1, 2) in c2
    assert El(1, 3) not in c2


def test_Diagram_name():
    assert assoc.name == 'assoc'


def test_Diagram_shape():
    binary = Shape.simplex(2).dual()
    assoc_l = binary.to_inputs(0, binary)
    assoc_r = binary.to_inputs(1, binary)
    assert assoc.shape == assoc_l.atom(assoc_r)


def test_Diagram_ambient():
    assert a.ambient == Mon
    assert c1.ambient == RP3


def test_Diagram_mapping():
    assert c2.mapping == [
            ['c0', 'c0', 'c0'],
            ['c0', 'c1', 'c1'],
            ['c2']]


def test_Diagram_layers():
    diagram = u.paste(pt.unit()).paste(
            a.paste(u))
    assert diagram.layers == [
            u.paste(pt.unit()), a.paste(u)]


def test_Diagram_rewrite_steps():
    diagram = u.paste(pt.unit()).paste(
            a.paste(u))
    assert diagram.rewrite_steps == [
            pt.unit().paste(pt.unit()),
            a.paste(pt.unit()),
            a.paste(a)]


def test_Diagram_dim():
    assert c0.dim == 0
    assert c1.dim == 1
    assert c2.dim == 2
    assert c3.dim == 3


def test_Diagram_isdegenerate():
    assert not Diagram(C).isdegenerate
    assert a.lunitor('-').isdegenerate
    assert not c2.isdegenerate


def test_Diagram_isround():
    assert not m.paste(m, 0).isround
    assert m.paste(m, 0).paste(m).isround


def test_Diagram_iscell():
    assert m.iscell
    assert a.unit().iscell
    assert not m.paste(a).iscell


def test_Diagram_isinvertiblecell():
    assert a.lunitor('-').isinvertiblecell
    assert assoc.isinvertiblecell
    assert not m.isinvertiblecell


def test_Diagram_hascomposite():
    assert f.paste(g).hascomposite
    assert m.hascomposite
    assert not a.paste(a).hascomposite


def test_Diagram_rename():
    fpasteg = f.paste(g)
    assert str(fpasteg) == '(f) #0 (g)'
    fpasteg.rename('f #0 g')
    assert str(fpasteg) == 'f #0 g'


def test_Diagram_paste():
    with raises(ValueError) as err:
        f.paste(a)
    assert str(err.value) == utils.value_err(
            a, 'not the same ambient DiagSet')

    assert f.paste(g).input == f.input
    assert f.paste(g).output == g.output

    assert f.paste(y, 0) == f
    assert x.paste(f, 0) == f

    with raises(ValueError) as err:
        f.paste(f)
    assert str(err.value) == utils.value_err(
            f, 'boundary does not match boundary of {}'.format(
                repr(f)))


def test_Diagram_to_outputs():
    with raises(ValueError) as err:
        f.to_outputs(1, a)
    assert str(err.value) == utils.value_err(
            a, 'not the same ambient DiagSet')

    with raises(ValueError) as err:
        c2.to_outputs(1, c2)
    assert str(err.value) == utils.value_err(
            c2, 'boundary does not match boundary of {}'.format(
                repr(c2)))

    assert c2.to_outputs(1, c1.unit()).output == c2.output


def test_Diagram_to_inputs():
    with raises(ValueError) as err:
        f.to_inputs(0, a)
    assert str(err.value) == utils.value_err(
            a, 'not the same ambient DiagSet')

    pt2unit = pt.unit().unit()
    with raises(ValueError) as err:
        m.to_inputs(0, pt2unit)
    assert str(err.value) == utils.value_err(
            m,
            'boundary does not match boundary of {}'.format(
                repr(pt2unit)))

    assert m.to_inputs(1, m).input == a.paste(a).paste(a)


def test_Diagram_pullback():
    arrow = Shape.arrow()
    connection = arrow.cube_connection(0, '-')

    fconn = f.pullback(connection)
    assert fconn.shape == connection.source

    with raises(ValueError) as err:
        m.pullback(connection)
    assert str(err.value) == utils.value_err(
            connection, 'target does not match diagram shape')


def test_Diagram_boundary():
    assert m.input == a.paste(a)
    assert m.output == a
    assert m.boundary('-', 0) == pt
    assert m.boundary('+', 0) == pt


def test_Diagram_unit():
    diagram = m.to_inputs(0, u)
    unit = diagram.unit()
    assert unit.shape == diagram.shape.inflate().source
    assert unit.input == diagram
    assert unit.output == diagram


def test_Diagram_lunitor():
    assert a.lunitor('-').input == pt.unit().paste(a)
    assert a.lunitor('+') == a.lunitor('-').inverse
    assert m.lunitor('-', 0).input == m.to_inputs(0, a.unit())


def test_Diagram_runitor():
    assert a.runitor('+').output == a.paste(pt.unit())
    assert a.runitor('+') == a.runitor('-').inverse
    assert c2.runitor('-', 2).input == c2.to_outputs(2, c1.unit())


def test_Diagram_inverse():
    assert assoc.inverse == assinv
    assert a.unit().inverse == a.unit()
    with raises(ValueError) as err:
        m.inverse
    assert str(err.value) == utils.value_err(
            m, 'not an invertible cell')


def test_Diagram_rinvertor():
    assert assoc.rinvertor == assinv.linvertor
    assert assoc.rinvertor == rinv
    assert a.unit().rinvertor == a.unit().lunitor('-')
    with raises(ValueError) as err:
        m.rinvertor
    assert str(err.value) == utils.value_err(
            m, 'not an invertible cell')


def test_Diagram_linvertor():
    assert assoc.linvertor == linv
    assert a.unit().linvertor == a.unit().lunitor('-')
    with raises(ValueError) as err:
        m.linvertor
    assert str(err.value) == utils.value_err(
            m, 'not an invertible cell')


def test_Diagram_composite():
    assert f.paste(g).composite == h
    assert m.composite == m

    aa = a.paste(a)
    with raises(ValueError) as err:
        aa.composite
    assert str(err.value) == utils.value_err(
        aa, 'does not have a composite')


def test_Diagram_compositor():
    assert f.paste(g).compositor == c_fg
    assert m.compositor == m.unit()

    aa = a.paste(a)
    with raises(ValueError) as err:
        aa.compositor
    assert str(err.value) == utils.value_err(
        aa, 'does not have a compositor')


def test_Diagram_yoneda():
    arrow = Shape.arrow()
    connection = arrow.cube_connection(0, '-')
    yoneda_conn = Diagram.yoneda(connection)

    assert yoneda_conn.ambient.generators == \
        DiagSet.yoneda(arrow).generators
    assert yoneda_conn.shape == connection.source
    assert yoneda_conn.mapping == connection.mapping


def test_Diagram_with_layers():
    layer1 = m.paste(pt.unit())
    layer2 = a.paste(u)
    layer3 = m

    diagram = Diagram.with_layers(layer1, layer2, layer3)
    assert diagram.layers == [layer1, layer2, layer3]
    assert diagram == layer1.paste(
            layer2.paste(layer3))


""" Tests for Diagram subclasses """


def test_SimplexDiagram():
    assert isinstance(c3.simplex_face(2), diagrams.SimplexDiagram)
    assert isinstance(
            c1.simplex_degeneracy(1), diagrams.SimplexDiagram)

    assert c1.simplex_degeneracy(1).simplex_degeneracy(2) == \
        c1.simplex_degeneracy(1).simplex_degeneracy(1)
    assert c2.simplex_degeneracy(2).simplex_face(2) == c2
    assert c2.simplex_degeneracy(2).simplex_face(0) == \
        c2.simplex_face(0).simplex_degeneracy(1)


def test_CubeDiagram():
    assert isinstance(s.cube_face(1, '+'), diagrams.CubeDiagram)
    assert isinstance(e1.cube_degeneracy(1), diagrams.CubeDiagram)
    assert isinstance(e0.cube_connection(0, '+'), diagrams.CubeDiagram)

    assert e0.cube_degeneracy(1).cube_degeneracy(2) == \
        e0.cube_degeneracy(1).cube_degeneracy(1)
    assert e1.cube_connection(0, '-').cube_connection(1, '-') == \
        e1.cube_connection(0, '-').cube_connection(0, '-')
