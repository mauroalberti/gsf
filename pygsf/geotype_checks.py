# -*- coding: utf-8 -*-


from .ptbaxes import *


def is_gvect(obj) -> bool:

    return isinstance(obj, GVect) and not isinstance(obj, GAxis)


def is_gaxis(obj) -> bool:

    return isinstance(obj, GAxis)


def is_gplane(obj) -> bool:

    return isinstance(obj, GPlane)


def is_slickln(obj) -> bool:

    return isinstance(obj, Slick)


def is_faultslck(obj) -> bool:

    return isinstance(obj, GFault)


def is_ptbaxes(obj) -> bool:

    return isinstance(obj, PTBAxes)


def is_upward(obj: GVect) -> bool:

    return obj.isUpward


def is_not_upward(obj: GVect) -> bool:

    return not obj.isUpward



