from rewal.shapes import Shape


def test_Shape():
    assert Shape() == Shape()
    assert Shape().size == [1]
