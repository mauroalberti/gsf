# -*- coding: utf-8 -*-


class AnaliticSurfaceCalcException(Exception):
    """
    Exception for Analytical Surface calculation.
    """

    pass


class SubparallelLineationException(Exception):
    """
    Exception for subparallel GAxis/GVect instances.
    """

    pass


class SlickelineTypeException(Exception):
    """
    Exception for slickenline type.
    """

    pass


class SlickelineSenseException(Exception):
    """
    Exception for slickenline movement sense.
    """

    pass


if __name__ == "__main__":

    import doctest
    import numtest  # external module, used in doctest float checks
    doctest.testmod()

