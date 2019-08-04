# -*- coding: utf-8 -*-

import numbers


min_epsg_crs_code = 2000  # checked 2019-06-14 in EPSG database


class Crs(object):
    """
    CRS class.
    Currently it is in a basic form,
    just managing simple comparisons and validity checks.

    """

    def __init__(self, epsg_cd: numbers.Integral = -1):

        self._epsg = int(epsg_cd)

    def epsg(self) -> numbers.Integral:

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


def check_crs(
    el1,
    el2
) -> bool:
    """
    Check whether two spatial elements have the same crs.

    :param el1: first spatial element.
    :param el2: second spatial element.
    :param err_msg: the error msg.
    :type err_msg: str.
    :return: whether two spatial elements have the same crs.
    :rtype: bool.
    """

    if el1.crs() != el2.crs():
        raise Exception("First {} has {} EPSG code while second {} has {}".format(
            type(el1),
            el1.epsg(),
            type(el2),
            el2.epsg()
        )
    )
