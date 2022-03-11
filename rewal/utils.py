"""
Utility functions for rewal.
"""


def typecheck(x, constraint, *more_constraints):
    """
    Type & constraint checking function.
    """
    if not isinstance(x, constraint['type']):
        raise TypeError(type_err(constraint['type'], x))
    if 'st' in constraint:
        if not constraint['st'](x):
            raise ValueError(value_err(x, constraint['why']))

    if more_constraints and hasattr(x, '__iter__'):
        for y in x:
            if isinstance(x, dict):
                typecheck(x[y], *more_constraints)
            else:
                typecheck(y, *more_constraints)


def type_err(expected, got):
    """ Type error. """
    return "Expected {}.{}, got {} of type {}.{} instead.".format(
           expected.__module__, expected.__name__,
           repr(got), type(got).__module__, type(got).__name__)


def value_err(got, why):
    """ Value error. """
    return "{} is not a valid value ({}).".format(repr(got), why)


def make_sign(key):
    """ Used to turn various expressions into '-' or '+' """
    if isinstance(key, str):
        if key in ['-', 'i', 'in', 'input', 's', 'source']:
            return '-'
        if key in ['+', 'o', 'out', 'output', 't', 'target']:
            return '+'
    if isinstance(key, int):
        if key == 0:
            return '-'
        if key == 1:
            return '+'
    raise KeyError(str(key))


def flip(sign):
    """ Flips the sign. """
    flipped = '-' if sign == '+' else '+'
    return flipped
