

def string2dict(str, valsep=",", kvsep="="):
    """
    Creates a dictionary from a string.

    :param str: string to convert into dictionary
    :param valsep: separator between key-value pairs
    :param kvsep: separator between key and value
    :return: a dictionary

    Examples:
      >>> string2dict("m=s, c=blue")
      {'c': 'blue', 'm': 's'}
    """

    vals = str.split(valsep)
    kv_vals = map(lambda kvstr: kvstr.strip().split(kvsep), vals)
    return dict(kv_vals)


