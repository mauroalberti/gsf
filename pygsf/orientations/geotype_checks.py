# -*- coding: utf-8 -*-


from pygsf.geology.structural.ptbaxes import *


def isOrien(obj) -> bool:

    return isinstance(obj, Orien) and not isinstance(obj, Axis)


def isAxis(obj) -> bool:

    return isinstance(obj, Axis)


def isPPlane(obj) -> bool:

    return isinstance(obj, PPlane)


def isSlickln(obj) -> bool:

    return isinstance(obj, Slick)


def isFaultSlck(obj) -> bool:

    return isinstance(obj, GFault)


def isPTBAxes(obj) -> bool:

    return isinstance(obj, PTBAxes)


def is_upward(obj: Orien) -> bool:

    return obj.isUpward


def is_not_upward(obj: Orien) -> bool:

    return not obj.isUpward



