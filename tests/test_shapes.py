import numpy as np
from pytest import raises

from rewal import utils
from rewal.ogposet import El, OgPoset, GrSet, GrSubset, Closed, OgMap


""" Tests for El """


def test_El_init():
    assert El(2, 3) == El(2, 3)
    assert El(2, 3) != El(1, 3)
    assert El(2, 3) != El(2, 2)

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
    assert str(err.value) == utils.value_err(
            -4, 'shifted position must be non-negative')


""" Tests for OgPoset """


def test_OgPoset_init():
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

    assert OgPoset(test_face, test_coface) == OgPoset(test_face, test_coface)

    test_face[0][0]['-'] = {0}
    with raises(ValueError) as err:
        OgPoset(test_face, test_coface)
    assert str(err.value) == utils.value_err(0, 'out of bounds')

    test_face[0][0]['-'] = set()
    test_face[1][0]['-'] = {2}
    with raises(ValueError) as err:
        OgPoset(test_face, test_coface)
    assert str(err.value) == utils.value_err(2, 'out of bounds')

    test_face[1][0]['-'] = {1}
    with raises(ValueError) as err:
        OgPoset(test_face, test_coface)
    assert str(err.value) == utils.value_err(
            test_face, 'input and output faces of El(1, 0) are not disjoint')

    test_face[1][0]['-'] = set()
    test_face[1][0]['+'] = set()
    with raises(ValueError) as err:
        OgPoset(test_face, test_coface)
    assert str(err.value) == utils.value_err(
            test_face, 'El(1, 0) must have at least one face')

    test_face[1][0]['-'] = {1}
    test_face[1][0]['+'] = {0}
    with raises(ValueError) as err:
        OgPoset(test_face, test_coface)
    assert str(err.value) == 'Face and coface data do not match.'

    test_face[1][0]['-'] = {'x'}
    with raises(TypeError) as err:
        OgPoset(test_face, test_coface)
    assert str(err.value) == utils.type_err(int, 'x')

    test_face[1][0]['-'] = 0
    with raises(TypeError) as err:
        OgPoset(test_face, test_coface)
    assert str(err.value) == utils.type_err(set, 0)

    test_face[1][0] = {'k': {0}}
    with raises(ValueError) as err:
        OgPoset(test_face, test_coface)
    assert str(err.value) == utils.value_err(
            {'k': {0}},
            "expecting dict with keys '-', '+'")

    test_face[1][0] = {0}
    with raises(TypeError) as err:
        OgPoset(test_face, test_coface)
    assert str(err.value) == utils.type_err(dict, {0})

    test_face[1] = []
    with raises(ValueError) as err:
        OgPoset(test_face, test_coface)
    assert str(err.value) == utils.value_err(
            [], 'expecting non-empty list')

    test_face[1] = {0}
    with raises(TypeError) as err:
        OgPoset(test_face, test_coface)
    assert str(err.value) == utils.type_err(list, {0})

    test_face = {0}
    with raises(TypeError) as err:
        OgPoset(test_face, test_coface)
    assert str(err.value) == utils.type_err(list, {0})


""" Various example objects here """

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

whisker = OgPoset(whisker_face, whisker_coface)

interval_face = [
        [
            {'-': set(), '+': set()},
            {'-': set(), '+': set()},
        ], [
            {'-': {0}, '+': {1}},
        ]]

interval = OgPoset.from_face_data(interval_face)

test_grset = GrSet(El(0, 2), El(2, 0), El(0, 5))

test_grsubset = GrSubset(GrSet(El(0, 2), El(2, 0)), whisker)
test_closed = test_grsubset.closure()

interval_grsubset = GrSubset(GrSet(El(0, 1)), interval)

whisker_all = whisker.all_elements

test_injection = OgMap(interval, whisker, [
    [El(0, 1), El(0, 2)],
    [El(1, 2)]])

test_collapse = OgMap(whisker, interval, [
    [El(0, 0), El(0, 1), El(0, 1)],
    [El(1, 0), El(1, 0), El(0, 1)],
    [El(1, 0)]])

test_composite = OgMap(interval, interval, [
    [El(0, 1), El(0, 1)],
    [El(0, 1)]])


def test_OgPoset_str():
    assert str(whisker) == 'OgPoset with [3, 3, 1] elements'


def test_OgPoset_getitem():
    assert whisker[2] == GrSubset(GrSet(El(2, 0)), whisker)


def test_OgPoset_contains():
    assert El(1, 2) in whisker
    assert El(1, 3) not in whisker
    assert El(3, 0) not in whisker


def test_OgPoset_len():
    assert len(whisker) == 7
    assert len(interval) == 3


def test_OgPoset_size():
    assert whisker.size == [3, 3, 1]


def test_OgPoset_dim():
    assert whisker.dim == 2


def test_OgPoset_as_chain():
    chain = [
            np.array([[-1, -1, 0], [1, 1, -1], [0, 0, 1]]),
            np.array([[-1], [1], [0]])
            ]
    test_chain = whisker.as_chain
    assert (test_chain[0] == chain[0]).all() and \
        (test_chain[1] == chain[1]).all()


def test_OgPoset_all_elements():
    assert whisker.all_elements == Closed(
            GrSet(El(0, 0), El(0, 1), El(0, 2),
                  El(1, 0), El(1, 1), El(1, 2),
                  El(2, 0)),
            whisker)


def test_OgPoset_faces():
    assert whisker.faces(El(2, 0), '-') == GrSet(El(1, 0))
    assert whisker.faces(El(2, 0), '+') == GrSet(El(1, 1))
    assert whisker.faces(El(0, 0), 0) == GrSet()
    assert whisker.faces(El(1, 2)) == GrSet(El(0, 1), El(0, 2))


def test_OgPoset_from_face_data():
    assert whisker == OgPoset.from_face_data(whisker_face)


""" Tests for GrSet """


def test_GrSet_init():
    assert GrSet(El(0, 2), El(1, 4), El(0, 3)) \
            == GrSet(El(0, 2), El(1, 4), El(0, 3))

    with raises(TypeError) as err:
        GrSet((0, 2))
    assert str(err.value) == utils.type_err(El, (0, 2))


def test_GrSet_str():
    assert str(test_grset) == 'GrSet[El(0, 2), El(0, 5), El(2, 0)]'


def test_GrSet_contains():
    assert El(0, 5) in test_grset
    assert El(0, 4) not in test_grset
    assert El(1, 0) not in test_grset
    assert 'x' not in test_grset


def test_GrSet_len():
    assert len(test_grset) == 3


def test_GrSet_iter():
    assert test_grset == GrSet(*test_grset)


def test_GrSet_getitem():
    assert test_grset[0] == GrSet(El(0, 2), El(0, 5))
    assert test_grset[1] == GrSet()

    assert test_grset[:3] == test_grset
    assert test_grset[1:] == test_grset[2]
    assert test_grset[:] == test_grset

    with raises(KeyError) as err:
        test_grset['x']
    assert str(err.value) == "'x'"


def test_GrSet_grades():
    assert test_grset.grades == [0, 2]


def test_GrSet_dim():
    assert test_grset.dim == 2
    assert GrSet().dim == -1


def test_GrSet_as_set():
    assert test_grset.as_set == \
            {El(0, 2), El(2, 0), El(0, 5)}


def test_GrSet_as_list():
    assert test_grset.as_list == \
            [El(0, 2), El(0, 5), El(2, 0)]


def test_GrSet_add():
    test_grset.add(El(3, 6))
    assert test_grset == GrSet(El(0, 2), El(3, 6), El(2, 0), El(0, 5))
    test_grset.add(El(0, 5))
    assert test_grset == GrSet(El(0, 2), El(3, 6), El(2, 0), El(0, 5))

    with raises(TypeError) as err:
        test_grset.add((3, 5))
    assert str(err.value) == utils.type_err(El, (3, 5))


def test_GrSet_remove():
    test_grset.remove(El(3, 6))
    assert test_grset == GrSet(El(0, 2), El(2, 0), El(0, 5))

    with raises(ValueError) as err:
        test_grset.remove(El(3, 6))
    assert str(err.value) == utils.value_err(
            El(3, 6), 'not in graded set')


def test_GrSet_union():
    assert test_grset.union(GrSet(El(1, 3), El(2, 0))) == \
            GrSet(El(0, 2), El(0, 5), El(1, 3), El(2, 0))


def test_GrSet_intersection():
    assert test_grset.intersection(GrSet(El(1, 3), El(2, 0))) == \
            GrSet(El(2, 0))
    assert test_grset.intersection(GrSet(El(1, 3))) == GrSet()


def test_GrSet_issubset():
    assert test_grset.issubset(test_grset)
    assert GrSet(El(0, 2), El(2, 0)).issubset(test_grset)
    assert not GrSet(El(0, 3)).issubset(test_grset)


def test_GrSet_isdisjoint():
    assert test_grset.isdisjoint(GrSet(El(1, 3), El(0, 4)))
    assert not test_grset.isdisjoint(GrSet(El(0, 2), El(1, 3)))


""" Tests for GrSubset """


def test_GrSubset():
    assert GrSubset(GrSet(El(0, 2), El(2, 0)), whisker) == \
        GrSubset(GrSet(El(0, 2), El(2, 0)), whisker)
    assert GrSubset(GrSet(El(0, 1)), interval) != \
        GrSubset(GrSet(El(0, 1)), whisker)
    assert GrSubset(GrSet(El(0, 1)), whisker) != \
        GrSubset(GrSet(El(0, 2)), whisker)


def test_GrSubset_init():
    with raises(ValueError) as err:
        GrSubset(test_grset, whisker)
    assert str(err.value) == utils.value_err(
            test_grset, 'does not define a subset')


def test_GrSubset_str():
    assert str(GrSubset(GrSet(El(0, 1)), interval)) == \
            'GrSubset with 1 elements in OgPoset with [2, 1] elements'


def test_GrSubset_contains():
    assert El(0, 2) in test_grsubset
    assert El(1, 1) not in test_grsubset


def test_GrSubset_getitem():
    assert test_grsubset[0] == GrSubset(GrSet(El(0, 2)), whisker)
    assert test_grsubset[1:] == GrSubset(GrSet(El(2, 0)), whisker)


def test_GrSubset_proj():
    assert test_grsubset.proj == GrSet(El(0, 2), El(2, 0))


def test_GrSubset_ambient():
    assert test_grsubset.ambient == whisker


def test_GrSubset_isclosed():
    assert not test_grsubset.isclosed
    assert whisker_all.isclosed


def test_GrSubset_union():
    assert test_grsubset.union(GrSubset(
        GrSet(El(0, 2), El(1, 1)), whisker), GrSubset(
        GrSet(El(2, 0), El(1, 2)), whisker)) == GrSubset(
        GrSet(El(0, 2), El(2, 0), El(1, 1), El(1, 2)), whisker)

    with raises(ValueError) as err:
        test_grsubset.union(interval_grsubset)
    assert str(err.value) == utils.value_err(
            interval_grsubset,
            'not a subset of the same OgPoset')

    assert not isinstance(
            test_grsubset.union(whisker.all_elements),
            Closed)
    assert isinstance(
            test_grsubset.closure().union(whisker.all_elements),
            Closed)


def test_GrSubset_intersection():
    assert test_grsubset.intersection(GrSubset(
        GrSet(El(0, 2), El(1, 1)), whisker)) == GrSubset(
        GrSet(El(0, 2)), whisker)

    assert not isinstance(
            test_grsubset.intersection(whisker.all_elements),
            Closed)
    assert isinstance(
            test_grsubset.closure().intersection(whisker.all_elements),
            Closed)


def test_GrSubset_closure():
    assert test_grsubset.closure() == Closed(
            GrSet(El(0, 0), El(0, 1), El(0, 2),
                  El(1, 0), El(1, 1), El(2, 0)),
            whisker)
    assert whisker_all.maximal().closure() == \
        whisker_all


def test_GrSubset_image():
    assert test_grsubset.image(test_collapse) == GrSubset(
            GrSet(El(0, 1), El(1, 0)), interval)

    with raises(ValueError) as err:
        test_grsubset.image(test_injection)
    assert str(err.value) == utils.value_err(
            test_injection, 'OgMap source does not match ambient OgPoset')

    assert isinstance(
            test_closed.image(test_collapse),
            Closed)
    assert not isinstance(
            test_grsubset.image(test_collapse),
            Closed)


""" Tests for Closed """


def test_Closed():
    assert Closed(
            GrSet(El(0, 0), El(0, 1), El(1, 1)),
            whisker) == GrSubset(
                    GrSet(El(1, 1)), whisker).closure()

    with raises(ValueError) as err:
        Closed(GrSet(El(1, 1)), whisker)
    assert str(err.value) == utils.value_err(
            GrSet(El(1, 1)), 'not a closed subset')


def test_Closed_as_map():
    assert Closed(
            GrSet(El(0, 1), El(0, 2), El(1, 2)),
            whisker).as_map == test_injection

    assert test_closed.as_map.source.image(
            test_closed.as_map) == test_closed

    assert interval.image(test_injection).as_map == test_injection


def test_Closed_maximal():
    assert whisker_all.maximal() == GrSubset(
            GrSet(El(2, 0), El(1, 2)), whisker)

    assert test_closed.maximal() == test_grsubset


def test_Closed_boundary():
    assert whisker_all.boundary('-', 0) == Closed(
            GrSet(El(0, 0)), whisker)
    assert whisker_all.boundary('+', 0) == Closed(
            GrSet(El(0, 2)), whisker)
    assert whisker_all.boundary('s', 1) == GrSubset(
            GrSet(El(1, 0), El(1, 2)), whisker).closure()
    assert whisker_all.boundary('+') == GrSubset(
            GrSet(El(1, 1), El(1, 2)), whisker).closure()
    assert whisker_all.boundary(0, 2) == whisker_all

    assert whisker_all.boundary() == whisker_all[:2].closed()
    assert whisker_all.boundary(None, 0) == Closed(
            GrSet(El(0, 0), El(0, 2)), whisker)

    assert whisker_all.boundary('-', 3) == whisker_all
    assert whisker_all.boundary('-', -1) == Closed(
            GrSet(), whisker)

    assert Closed(GrSet(El(0, 0)), whisker).boundary('-') == \
        Closed(GrSet(), whisker)
    assert Closed(GrSet(), whisker).boundary('-') == \
        Closed(GrSet(), whisker)


""" Tests for OgMap """


def test_OgMap_init():
    assert OgMap(whisker, interval) == OgMap(whisker, interval)

    # TODO: tests for well-formedness


def test_OgMap_getitem():
    assert test_injection[El(0, 0)] == El(0, 1)
    assert test_collapse[El(2, 0)] == El(1, 0)
    assert OgMap(interval, whisker)[El(0, 1)] is None

    with raises(ValueError) as err:
        test_injection[El(0, 2)]
    assert str(err.value) == utils.value_err(
            El(0, 2), 'not in source')


def test_OgMap_setitem():
    test_setitem = OgMap(interval, whisker)
    test_setitem[El(0, 0)] = El(0, 1)

    with raises(ValueError) as err:
        test_setitem[El(0, 1)] = El(0, 3)
    assert str(err.value) == utils.value_err(
            El(0, 3), 'not in target')

    with raises(ValueError) as err:
        test_setitem[El(0, 0)] = El(0, 2)
    assert str(err.value) == utils.value_err(
            El(0, 0), 'already defined on element')

    with raises(ValueError) as err:
        test_setitem[El(0, 1)] = El(1, 0)
    assert str(err.value) == utils.value_err(
            El(1, 0), 'exceeds dimension of El(0, 1)')

    with raises(ValueError) as err:
        test_setitem[El(1, 0)] = El(1, 2)
    assert str(err.value) == utils.value_err(
            El(1, 0), 'map undefined on El(0, 1) below El(1, 0)')

    test_setitem[El(0, 1)] = El(0, 2)

    with raises(ValueError) as err:
        test_setitem[El(1, 0)] = El(1, 1)
    assert str(err.value) == utils.value_err(
            El(1, 1),
            'assignment does not respect (-, 0)-boundary')

    test_setitem[El(1, 0)] = El(1, 2)

    assert test_setitem == test_injection


def test_OgMap_mapping():
    assert OgMap(interval, whisker).mapping == [[None, None], [None]]


def test_OgMap_istotal():
    assert test_injection.istotal
    assert not OgMap(interval, whisker).istotal


def test_OgMap_isinjective():
    assert test_injection.isinjective
    assert not test_collapse.isinjective


def test_OgMap_issurjective():
    assert not test_injection.issurjective
    assert test_collapse.issurjective


def test_OgMap_isdefined():
    assert test_injection.isdefined(El(0, 1))
    assert not OgMap(interval, whisker).isdefined(El(0, 1))
    assert not test_injection.isdefined(El(0, 2))


def test_OgMap_then():
    assert test_injection.then(test_collapse) == test_composite
    assert test_injection.then(OgMap(whisker, whisker)) == \
        OgMap(interval, whisker)
    assert test_injection.then(
            test_collapse, OgMap(interval, interval)) == \
        OgMap(interval, interval)

    with raises(ValueError) as err:
        test_injection.then(test_injection)
    assert str(err.value) == utils.value_err(
            test_injection, 'source does not match target of first map')
