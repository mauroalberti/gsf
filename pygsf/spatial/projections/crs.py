# -*- coding: utf-8 -*-


min_epsg_crs_code = 2000  # checked 2019-06-14 in EPSG database


class Crs(object):
    """
    CRS class.
    Currently it is in a basic form,
    just managing simple comparisons and validity checks.

    """

    def __init__(self, epsg_cd: int = -1):

        self._epsg = int(epsg_cd)

    def epsg(self) -> int:

        return self._epsg

    def valid(self):

        return self.epsg() >= min_epsg_crs_code

    def __repr__(self):

        return "EPSG:{}".format(self.epsg())

    def __eq__(self, another) -> bool:
        """
        Checks for equality between Crs instances.
        Currently it considers equal two Crs instances when they have the
        same EPSG code, even an invalid one (i.e., -1).

        :param another: the Crs instance to compare with.
        :type another: Crs.
        :return: whether the input Crs instance is equal to the current one.
        :rtype: bool.
        :raise: Exception.
        """

        if not (isinstance(another, Crs)):
            raise Exception("Input instance should be Crs but is {}".format(type(another)))

        return self.epsg() == another.epsg()


