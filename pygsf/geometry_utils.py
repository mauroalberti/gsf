

def mod360(val: float) -> float:
    """
    Return the module 360 of the originila value.
    
    :param val: value
    :type val: float
    :return: value divided by mod 360
    :rtype: float

    Examples:
    >>> mod360(10)
    10.0
    >>> mod360(360)
    0.0
    >>> mod360(-50)
    310.0
    >>> mod360(400)
    40.0

    """

    return val % 360.0


def opposite_trend(tr: float) -> float:
    """
    Calculate the trend opposite to the original one.
    
    :return: the opposite trend.

    Examples:
    >>> opposite_trend(0)
    180.0
    >>> opposite_trend(45)
    225.0
    >>> opposite_trend(90)
    270.0
    >>> opposite_trend(180)
    0.0
    >>> opposite_trend(270)
    90.0
    """
    
    return mod360(tr + 180.0)


def plng2colatTop(plunge: float) -> float:
    """
    Calculates the colatitude angle from the top.

    :param plunge: an angle from -90° (upward-pointing) to 90° (downward-pointing)
    :type plunge: float
    :return: the colatitude angle
    :rtype: float

    Examples:
      >>> plng2colatTop(90)
      180.0
      >>> plng2colatTop(45)
      135.0
      >>> plng2colatTop(0)
      90.0
      >>> plng2colatTop(-45)
      45.0
      >>> plng2colatTop(-90)
      0.0
    """

    return 90.0 + plunge


def plng2colatBottom(plunge: float) -> float:
    """
    Calculates the colatitude angle from the bottom.

    :param plunge: an angle from -90° (upward-pointing) to 90° (downward-pointing)
    :type plunge: float
    :return: the colatitude angle
    :rtype: float

    Examples:
      >>> plng2colatBottom(90)
      0.0
      >>> plng2colatBottom(45)
      45.0
      >>> plng2colatBottom(0)
      90.0
      >>> plng2colatBottom(-45)
      135.0
      >>> plng2colatBottom(-90)
      180.0
    """

    return 90.0 - plunge



if __name__ == "__main__":

    import doctest
    doctest.testmod()
