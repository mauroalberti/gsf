
from ..ptbaxes import *


def isDirect(obj) -> bool:

    return isinstance(obj, Direct) and not isinstance(obj, Axis)


def isAxis(obj) -> bool:

    return isinstance(obj, Axis)


def isPlane(obj) -> bool:

    return isinstance(obj, Plane)


def isSlick(obj) -> bool:

    return isinstance(obj, Slick)


def isFault(obj) -> bool:

    return isinstance(obj, Fault)


def isPTBAxes(obj) -> bool:

    return isinstance(obj, PTBAxes)





