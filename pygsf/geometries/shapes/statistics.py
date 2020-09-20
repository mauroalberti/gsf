
from typing import Tuple, Union, List

from functools import singledispatch
import numpy as np

from pygsf.mathematics.arrays import xyzSvd
from pygsf.mathematics.vectors import Vect
from pygsf.orientations.orientations import Direct
from pygsf.geometries.shapes.space2d import *
from pygsf.geometries.shapes.space3d import CPlane
from pygsf.geolocated.geoshapes import GeoPoints


@singledispatch
def mean(
        shapes: List[Shape2D]
) -> Shape2D:

    return None


@mean.register(List[Point])
def mean(
        shapes: List[Point]
) -> Point:
    """Mean points center"""

    return Point(
        x=np.mean(shapes.xs()),
        y=np.mean(shapes.ys())
    )


def try_derive_bestfitplane(
    points: GeoPoints
) -> Tuple[bool, Union[str, CPlane]]:

    print(points.xs)
    print(points.ys)
    print(points.zs)

    npaXyz = points.asXyzArray()

    print(points.asXyzArray())

    xyz_mean = np.mean(npaXyz, axis=0)

    svd = xyzSvd(npaXyz - xyz_mean)

    if svd['result'] is None:
        return False, "Unable to calculate result"

    _, _, eigenvectors = svd['result']

    lowest_eigenvector = eigenvectors[-1, : ]  # Solution is last row

    normal = lowest_eigenvector[: 3 ] / np.linalg.norm(lowest_eigenvector[: 3 ])
    normal_vector = Vect(normal[0], normal[1], normal[2])
    normal_direct = Direct.fromVect(normal_vector)

    return True, normal_direct.normal_plane()