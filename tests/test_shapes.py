from pytest import raises

from rewalt import utils, shapes
from rewalt.ogposets import (El, OgPoset, OgMap, OgMapPair)
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


def test_Shape_to_inputs():
    arrow = Shape.arrow()
    twothree = arrow.paste(arrow).atom(arrow.paste(arrow).paste(arrow))
    threetwo = twothree.dual()

    with raises(ValueError) as err:
        threetwo.to_inputs([0, 2], threetwo)
    assert str(err.value) == utils.value_err(
            [0, 2], 'cannot paste to these inputs')

    with raises(ValueError) as err:
        threetwo.to_inputs([0, 1], twothree)
    assert str(err.value) == utils.value_err(
            [0, 1], 'does not match output boundary of {}'.format(
                repr(twothree)))

    pasted = threetwo.to_inputs([0, 1], threetwo)
    assert pasted.size == [7, 8, 2]


def test_Shape_suspend():
    globe2 = Shape.globe(2)
    arrow = Shape.arrow()

    whisker_l = arrow.paste(globe2)

    assert whisker_l.suspend().size == [2] + whisker_l.size
    assert arrow.suspend(2) == globe2.suspend()

    assert isinstance(arrow.paste(arrow).suspend(), shapes.GlobeString)
    assert isinstance(whisker_l.suspend(), shapes.Theta)


def test_Shape_gray():
    globe2 = Shape.globe(2)
    arrow = Shape.arrow()

    assert (globe2 * arrow).size == [4, 6, 4, 1]
    assert Shape.gray() == Shape.point()
    assert (arrow * arrow) * arrow == arrow * (arrow * arrow)

    assert isinstance(arrow * arrow, shapes.Cube)


def test_Shape_join():
    point = Shape.point()
    arrow = Shape.arrow()

    assert arrow >> point == point << point << point
    assert Shape.join() == Shape.empty()
    assert (point >> point) >> point == point >> (point >> point)

    assert isinstance(arrow >> point, shapes.Simplex)


def test_Shape_dual():
    arrow = Shape.arrow()
    simplex = Shape.simplex(2)
    binary = arrow.paste(arrow).atom(arrow)
    assert binary == simplex.dual()

    assoc_l = binary.to_inputs(0, binary)
    assoc_r = binary.to_inputs(1, binary)
    assert assoc_r == assoc_l.dual(1)

    globe = Shape.globe(2)
    assert isinstance(arrow.paste(globe).dual(), shapes.Theta)


def test_Shape_merge():
    arrow = Shape.arrow()
    binary = arrow.paste(arrow).atom(arrow)
    ternary = arrow.paste(arrow).paste(arrow).atom(arrow)
    assoc_l = binary.to_inputs(0, binary)

    assert assoc_l.merge() == ternary
    assert isinstance(assoc_l.merge(), shapes.Opetope)


def test_Shape_empty():
    empty = Shape.empty()
    assert len(empty) == 0


def test_Shape_point():
    point = Shape.point()
    assert len(point) == 1


def test_Shape_arrow():
    arrow = Shape.arrow()
    assert arrow.size == [2, 1]


def test_Shape_simplex():
    assert len(Shape.simplex()) == 0
    arrow = Shape.simplex(1)
    assert arrow == Shape.arrow()
    assert arrow >> arrow == Shape.simplex(3)


def test_Shape_cube():
    assert len(Shape.cube()) == 1

    arrow = Shape.cube(1)
    assert arrow == Shape.arrow()
    assert arrow*arrow == Shape.cube(2)


def test_Shape_globe():
    assert len(Shape.globe()) == 1

    arrow = Shape.globe(1)
    assert arrow == Shape.arrow()
    assert arrow.suspend() == Shape.globe(2)


def test_Shape_theta():
    assert Shape.theta() == Shape.globe(0)
    assert Shape.theta(Shape.theta()) == Shape.globe(1)
    assert Shape.theta(Shape.theta(Shape.theta())) == Shape.globe(2)

    point = Shape.theta()
    arrow = Shape.arrow()
    assert Shape.theta(point, point) == arrow.paste(arrow)


def test_Shape_id():
    arrow = Shape.arrow()
    assert arrow.id().source == arrow
    assert arrow.id().target == arrow


def test_Shape_boundary():
    arrow = Shape.arrow()
    binary = arrow.paste(arrow).atom(arrow)

    assert binary.boundary('-').source == arrow.paste(arrow)
    assert binary.boundary('+').source == arrow
    assert binary.boundary('-', 0).source == Shape.point()
    assert binary.boundary('-').target == binary
    assert not isinstance(binary.boundary().source, Shape)

    assoc_l = binary.to_inputs(0, binary)
    assert isinstance(assoc_l.boundary('-').source, shapes.OpetopeTree)
    assert isinstance(assoc_l.boundary('+').source, shapes.Arrow)


def test_Shape_atom_inclusion():
    arrow = Shape.arrow()
    globe = Shape.globe(2)
    whisker_l = arrow.paste(globe)
    assert whisker_l.atom_inclusion(El(2, 0)).source == globe

    assert isinstance(
            whisker_l.atom_inclusion(El(2, 0)).source,
            shapes.Globe)
    binary = arrow.paste(arrow).atom(arrow)
    assoc_l = binary.to_inputs(0, binary)
    assert isinstance(
            assoc_l.atom_inclusion(El(2, 0)).source,
            shapes.Opetope)

    simplex = Shape.simplex(3)
    assert isinstance(
            simplex.atom_inclusion(El(2, 0)).source,
            shapes.Simplex)

    cube = Shape.cube(3)
    assert isinstance(
            cube.atom_inclusion(El(2, 0)).source,
            shapes.Cube)


def test_Shape_initial():
    point = Shape.point()
    empty = Shape.empty()
    assert point.initial() == empty.terminal()
    assert empty.initial() == empty.id()
    assert point.initial().istotal


def test_Shape_terminal():
    point = Shape.point()
    assert point.terminal() == point.id()


def test_Shape_inflate():
    simplex2 = Shape.simplex(2)
    assert simplex2.inflate().boundary('-').source == \
        simplex2
    whisker_l = Shape.arrow().paste(Shape.globe(2))
    assert whisker_l.inflate().boundary('+').source == \
        whisker_l


def test_Shape_all_layerings():
    globe = Shape.globe(2)
    chain = globe.paste(globe, 0)
    for n, x in enumerate(chain.all_layerings()):
        number = n+1
    assert number == 2


def test_Shape_generate_layering():
    arrow = Shape.arrow()
    globe = Shape.globe(2)
    chain = globe.paste(globe, 0)
    chain.generate_layering()
    assert chain.layers[0].source == arrow.paste(globe)
    assert chain.layers[1].source == globe.paste(arrow)
    chain.generate_layering()
    assert chain.layers[0].source == globe.paste(arrow)
    assert chain.layers[1].source == arrow.paste(globe)


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
    point = OgPoset.point()
    arrow = point >> point
    ogmap = OgMap(
            arrow, point,
            [[El(0, 0), El(0, 0)], [El(0, 0)]])
    with raises(TypeError) as err:
        ShapeMap(ogmap)
    assert str(err.value) == utils.type_err(Shape, arrow)

    arrow = Shape.arrow()
    ogmap = OgMap(
            arrow, point,
            [[El(0, 0), El(0, 0)], [El(0, 0)]])
    with raises(TypeError) as err:
        ShapeMap(ogmap)
    assert str(err.value) == utils.type_err(Shape, point)

    point = Shape.point()
    undefined = OgMap(point, arrow)
    with raises(ValueError) as err:
        ShapeMap(undefined)
    assert str(err.value) == utils.value_err(
            undefined,
            'a ShapeMap must be total')


def test_ShapeMap_then():
    point = Shape.point()
    arrow = Shape.arrow()
    terminal = arrow.terminal()
    first_inj = arrow.atom_inclusion(El(0, 0))
    assert isinstance(terminal.then(first_inj), ShapeMap)

    ogmap = OgMap(
            arrow, point,
            [[El(0, 0), El(0, 0)], [El(0, 0)]])
    assert not isinstance(ogmap.then(first_inj), ShapeMap)


def test_ShapeMap_layers():
    globe = Shape.globe(2)
    cospan = globe.paste(globe, cospan=True)
    twoglobes = cospan.target
    cospanterminal = cospan.then(twoglobes.terminal())
    assert twoglobes.terminal().layers == [
            cospanterminal.fst, cospanterminal.snd]


def test_ShapeMap_rewrite_steps():
    arrow = Shape.arrow()
    globe = Shape.globe(2)
    twoglobes = globe.paste(globe)
    assert twoglobes.terminal().rewrite_steps == [
            arrow.terminal(), arrow.terminal(), arrow.terminal()]


def test_ShapeMap_gray():
    point = Shape.point()
    arrow = Shape.arrow()
    terminal = arrow.terminal()
    first_inj = arrow.atom_inclusion(El(0, 0))
    assert ShapeMap.gray() == point.id()
    assert arrow.id() * arrow.id() == (arrow * arrow).id()
    assert (terminal * arrow.id()).then(
        first_inj * terminal) == terminal.then(first_inj) * terminal


def test_ShapeMap_join():
    empty = Shape.empty()
    arrow = Shape.arrow()
    terminal = arrow.terminal()
    first_inj = arrow.atom_inclusion(El(0, 0))
    assert ShapeMap.join() == empty.id()
    assert arrow.id() >> arrow.id() == (arrow >> arrow).id()
    assert (terminal >> arrow.id()).then(
        first_inj >> terminal) == terminal.then(first_inj) >> terminal


def test_ShapeMap_dual():
    arrow = Shape.arrow()
    degen0 = arrow.simplex_degeneracy(0)
    degen1 = arrow.simplex_degeneracy(1)
    assert degen0.op() == degen1

    cdegen0 = arrow.cube_degeneracy(0)
    cdegen1 = arrow.cube_degeneracy(1)
    assert cdegen0.co() == cdegen1
