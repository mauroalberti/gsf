
from ...types.utils import *
from .geometries import *


class Lines(list):
    """
    Collection of lines, inheriting from list.

    """

    def __init__(self, lines: List[Line]):

        check_type(lines, "Lines", List)
        for el in lines:
            check_type(el, "Line", Line)

        super(Lines, self).__init__(lines)


class MultiLines(list):
    """
    Collection of multilines, inheriting from list.

    """

    def __init__(self, multilines: List[MultiLine]):

        check_type(multilines, "MultiLines", List)
        for el in multilines:
            check_type(el, "MultiLine", MultiLine)

        super(MultiLines, self).__init__(multilines)

