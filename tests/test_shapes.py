from pytest import raises

from rewal import utils
from rewal.shapes import Shape


def test_Shape():
    assert Shape() == Shape()
    assert Shape().size == [1]
