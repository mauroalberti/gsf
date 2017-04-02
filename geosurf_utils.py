

def rhrstrk2dd(rhr_strk):
    """Converts RHR strike value to dip direction value.

    Example:
        >>> rhrstrk2dd(285.5)
        15.5
    """

    return (rhr_strk + 90.0) % 360.0


if __name__ == "__main__":
    import doctest
    doctest.testmod()