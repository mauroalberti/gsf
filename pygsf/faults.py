# -*- coding: utf-8 -*-

from .geometry import *
from .mathematics import are_close


class Slick(GVect):
    """
    Slickeline.
    It can be defined through a GVect instance, in which case it has a movement sense,
    or via a GAxis, when the movement sense is unknown or not sure.
    When the movement sense is known, the GVect instance indicates the displacement of the block that is:
    - for a horizontal or a dipping, non vertical fault: the upper block
    - for a vertical fault: the block individuated by the (formal) dip direction.
    """

    def __init__(self, trend: float, plunge: float, known=True):
        """"
        Slick constructor.
        The 'mov_lin' argument is a GVect or a GAxis instance. 
        Depending on that, the movement sense will be known
        when a GVect provided, and unknown/uncertain when a GAxis provided.

        Example:
          >>> Slick(90, 10)
          Slick(090.00, +10.00, True)
          #>>> Slick(90, 10, known=False)
          #Slick(090.00, +10.00, False)
        """

        """
        if not isinstance(trend, (int, float)):
            raise SlickInputTypeException("Trend must be a number")
        if not isinstance(plunge, (int, float)):
            raise SlickInputTypeException("Plunge must be a number")
        if not isinstance(known, bool):
            raise SlickInputTypeException("Known movement sense must be a boolean")
        """

        super().__init__(trend, plunge)

    def __repr__(self):

        known = not isinstance(self, GAxis)

        return "Slick({:06.2f}, {:+06.2f}, {})".format(self.tr, self.pl, known)


if __name__ == "__main__":

    import doctest
    doctest.testmod()
