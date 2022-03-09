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

    if more_constraints:
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
    return "The value {} is not valid ({}).".format(repr(got), why)
