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


def mksign(key):
    """ Used to turn various expressions into '-' or '+'. """
    if key in ['-', 0,
               'i', 'in', 'input',
               'd', 'dom', 'domain',
               's', 'src', 'source']:
        return '-'
    if key in ['+', 1,
               'c', 'cod', 'codomain',
               'o', 'out', 'output', 
               't', 'tgt', 'target']:
        return '+'
    raise KeyError(str(key))


def flip(sign):
    """ Flips the sign. """
    flipped = '-' if sign == '+' else '+'
    return flipped
