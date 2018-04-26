# -*- coding: utf-8 -*-


from pygsf.geology.ptbaxes import *


def isOrien(obj) -> bool:

    return isinstance(obj, Direct) and not isinstance(obj, Axis)


def isAxis(obj) -> bool:

    return isinstance(obj, Axis)


def isPPlane(obj) -> bool:

    return isinstance(obj, PPlane)


def isSlickln(obj) -> bool:

    return isinstance(obj, Slick)


def isFaultSlck(obj) -> bool:

    return isinstance(obj, Fault)


def isPTBAxes(obj) -> bool:

    return isinstance(obj, PTBAxes)


def is_upward(obj: Direct) -> bool:

    return obj.isUpward


def is_not_upward(obj: Direct) -> bool:

    return not obj.isUpward



