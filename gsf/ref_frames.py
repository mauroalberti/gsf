# -*- coding: utf-8 -*-

import numpy as np

from .mathematics import almost_zero, are_close


class RefFrame(object):

    def __init__(self, versor_x, versor_y, versor_z):

        assert versor_x.is_near_unit
        assert versor_y.is_near_unit
        assert versor_z.is_near_unit

        assert versor_x.is_suborthogonal(versor_y)
        assert versor_x.is_suborthogonal(versor_z)
        assert versor_y.is_suborthogonal(versor_z)

        self.i = versor_x
        self.j = versor_y
        self.k = versor_z

if __name__ == "__main__":

    import doctest
    doctest.testmod()
