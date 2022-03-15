from pytest import raises

from rewal import utils
from rewal.ogposets import (El, OgPoset, OgMap)
from rewal.shapes import (Shape, ShapeMap)


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
ogwhisker = OgPoset.from_face_data(whisker_face)

ogempty = OgPoset([], [])

ogpoint = OgPoset.from_face_data([
    [{'-': set(), '+': set()}]
    ])

interval_face = [
        [
            {'-': set(), '+': set()},
            {'-': set(), '+': set()},
        ], [
            {'-': {0}, '+': {1}},
        ]]
oginterval = OgPoset.from_face_data(interval_face)

oginjection = OgMap(oginterval, ogwhisker, [
    [El(0, 1), El(0, 2)],
    [El(1, 2)]])

ogcollapse = OgMap(ogwhisker, oginterval, [
    [El(0, 0), El(0, 1), El(0, 1)],
    [El(1, 0), El(1, 0), El(0, 1)],
    [El(1, 0)]])


""" Tests for Shape """


def test_Shape():
    assert Shape() == Shape()
    assert Shape().dim == -1


def test_Shape_point():
    assert Shape.point().size == [1]
    assert isinstance(Shape.point(), Shape)


def test_Shape_upgrade():
    assert not ogpoint == Shape.point()
    assert Shape._upgrade(ogpoint) == Shape.point()


""" Tests for ShapeMap """


def test_ShapeMap():
    assert ShapeMap(Shape(), Shape.point()) == \
            ShapeMap(Shape(), Shape.point())
    assert ShapeMap(Shape(), Shape.point()) != \
        OgMap(Shape(), Shape.point())

    with raises(TypeError) as err:
        ShapeMap(oginterval, ogwhisker)
    assert str(err.value) == utils.type_err(Shape, oginterval)


def test_ShapeMap_upgrade():
    assert not isinstance(ogcollapse, ShapeMap)
    assert isinstance(ShapeMap._upgrade(ogcollapse), ShapeMap)
