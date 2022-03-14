from rewal.shapes import Shape


def test_Shape():
    assert Shape() == Shape()
    assert Shape().dim == -1


def test_Shape_point():
    assert Shape.point().size == [1]
    assert isinstance(Shape.point(), Shape)
