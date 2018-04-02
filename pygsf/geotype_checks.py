# -*- coding: utf-8 -*-


from .ptbaxes import *


def is_gvect(obj) -> bool:

    return isinstance(obj, GVect) and not isinstance(obj, GAxis)


def is_gaxis(obj) -> bool:

    return isinstance(obj, GAxis)


def is_gplane(obj) -> bool:

    return isinstance(obj, GPlane)


def is_slickln(obj) -> bool:

    return isinstance(obj, Slickenline)


def is_faultslck(obj) -> bool:

    return isinstance(obj, FaultSlick)


def is_ptbaxes(obj) -> bool:

    return isinstance(obj, PTBAxes)


def is_upward(obj: GVect) -> bool:

    return obj.is_upward


def is_not_upward(obj: GVect) -> bool:

    return not obj.is_upward



