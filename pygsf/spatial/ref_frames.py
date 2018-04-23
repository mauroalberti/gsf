# -*- coding: utf-8 -*-


from .vector.vector import Vect


class RefFrame(object):

    def __init__(self, versor_x, versor_y):
        """
        Constructor.

        :param versor_x: Vect instance, representing x axis orientation
        :param versor_y: Vect instance, representing y axis orientation
        """

        if not (versor_x.isAlmostUnit and versor_y.isAlmostUnit):
            raise RefFrameInputException("Input vectors must be near unit")

        if not versor_x.isSubOrthogonal(versor_y):
            raise RefFrameInputException("Input vectors must be sub-orthogonal")

        self._x = versor_x
        self._y = versor_y

    @property
    def x(self):
        """
        Return the x as_vector,

        :return: Vect instance

        Example:
          >>> RefFrame(Vect(1,0,0), Vect(0,1,0)).x
          Vect(1.0000, 0.0000, 0.0000)
        """

        return self._x

    @property
    def y(self):
        """
        Return the y as_vector.

        :return: Vect instance

        Example:
          >>> RefFrame(Vect(1,0,0), Vect(0,1,0)).y
          Vect(0.0000, 1.0000, 0.0000)
        """

        return self._y

    @property
    def z(self):
        """
        Return the z as_vector.

        :return: Vect instance

        Example:
          >>> RefFrame(Vect(1,0,0), Vect(0,1,0)).z
          Vect(0.0000, 0.0000, 1.0000)
        """

        return self.x.vCross(self.y)


class RefFrameInputException(Exception):
    """
    Exception for RefFrame input
    """

    pass


if __name__ == "__main__":

    import doctest
    doctest.testmod()
