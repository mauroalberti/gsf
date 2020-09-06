
from collections import namedtuple

from pygsf.geology.faults import *
from pygsf.geometries.geom3d.shapes import Point
from pygsf.geometries.geom3d.abstract import CPlane


georef_att_flds = [
    'id',
    'posit',
    'attitude'
]

GeorefAttitude = namedtuple('GeorefAttitude', georef_att_flds)


class StructuralSet:

    def __init__(self,
            location: Point,
            stratifications: Optional[List[CPlane]] = None,
            foliations: Optional[List[CPlane]] = None,
            faults: Optional[List[Fault]] = None
    ):
        """
        Creates a structural set.

        :param location:
        :type location: Point.
        :param stratifications:
        :type stratifications: Optional[List[CPlane]].
        :param foliations:
        :type foliations: Optional[List[CPlane]].
        :param faults:
        :type faults: Optional[List[CPlane]].
        """

        if not isinstance(location, Point):
            raise Exception("Location should be Point but is {}".format(type(location)))

        checks = [
            (stratifications, "Stratification"),
            (foliations, "Foliations"),
            (faults, "Faults")
        ]

        for var, name in checks:
            if var:
                if not isinstance(var, List):
                    raise Exception("{} should be a List but is {}".format(name, type(var)))
                for el in var:
                    if not isinstance(el, CPlane):
                        raise Exception("{} should be CPlane but is {}".format(name, type(el)))

        self._location = location

        if not stratifications:
            self._strats = []
        else:
            self._strats = stratifications

        if not foliations:
            self._foliats = []
        else:
            self._foliats = foliations

        if not faults:
            self._faults = []
        else:
            self._faults = faults








