
from typing import Any, Optional, Tuple, Union
import numbers

import numpy as np

from pyproj.transformer import Transformer, TransformerGroup, AreaOfInterest


min_epsg_crs_code = 2000  # checked 2019-06-14 in EPSG database


class Crs(object):
    """
    CRS class.
    Currently it is in a basic form,
    just managing simple comparisons and validity checks.

    """

    def __init__(self, epsg_cd: numbers.Integral = -1):

        self._epsg = int(epsg_cd)

    @property
    def epsg_code(self) -> int:

        return self._epsg

    def valid(self):

        return self.epsg_code >= min_epsg_crs_code

    def __repr__(self):

        return "EPSG:{}".format(self.epsg_code)

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

        return self.epsg_code == another.epsg_code


def check_crs(
    template_element: Any,
    checked_element: Any
) -> None:
    """
    Check whether two spatial elements have the same georeferenced, raising an exception when not equal.
    The two elements should both implement the georeferenced property.

    :param template_element: first spatial element.
    :type template_element: Any
    :param checked_element: second spatial element.
    :type checked_element: Any
    :return: nothing
    :rtype: None
    :raise: Exception
    """

    if checked_element.crs != template_element.crs:
        raise Exception("checked {} instance has {} EPSG code but {} expected".format(
            type(checked_element).__name__,
            checked_element.epsg_code,
            template_element.epsg_code
        )
    )


def check_epsg(
    spatial_element: Any,
    epsg_code: numbers.Integral
) -> None:
    """
    Check whether a spatial element has a given EPSG code, raising an exception when not true.
    The spatial element should implement the epsg_code method.

    :param spatial_element: spatial element
    :type spatial_element: Any
    :param epsg_code: the EPSG code
    :type epsg_code: numbers.Integral
    :return: nothing
    :rtype: None
    :raise: Exception
    """

    if spatial_element.epsg_code != epsg_code:
        raise Exception("checked {} instance has {} EPSG code but {} expected".format(
            type(spatial_element).__name__,
            spatial_element.epsg_code,
            epsg_code
        )
    )


def project_xy(
    x: numbers.Real,
    y: numbers.Real,
    source_epsg_code: numbers.Integral,
    dest_epsg_code: numbers.Integral = 4326,
) -> Optional[Tuple[numbers.Real, numbers.Real]]:

    if source_epsg_code == -1:
        return None

    if dest_epsg_code == -1:
        return None

    transformer = Transformer.from_crs(
        f"epsg:{source_epsg_code}",
        f"epsg:{dest_epsg_code}"
    )

    x_prj, y_prj = transformer.transform(x, y)

    return x_prj, y_prj


def project_extent(
    x_min: numbers.Real,
    x_max: numbers.Real,
    y_min: numbers.Real,
    y_max: numbers.Real,
    source_epsg_code: numbers.Integral
) -> Optional[Tuple[numbers.Real, numbers.Real, numbers.Real, numbers.Real]]:

    result = project_xy(
        x=x_min,
        y=y_min,
        source_epsg_code=source_epsg_code
    )

    if result is None:
        return None

    lon_min, lat_min = result

    result = project_xy(
        x=x_max,
        y=y_max,
        source_epsg_code=source_epsg_code
    )

    if result is None:
        return None

    lon_max, lat_max = result

    return lon_min, lat_min, lon_max, lat_max


def try_project_xy_arrays(
    x_array: np.ndarray,
    y_array: np.ndarray,
    source_epsg_code: numbers.Integral,
    dest_epsg_code: numbers.Integral,
    area_of_interest: Tuple[numbers.Real, numbers.Real, numbers.Real, numbers.Real]
) -> Tuple[bool, Union[str, Tuple[np.ndarray, np.ndarray]]]:
    """
    WARNING: currently this method is experimental.
    To understand why axis swap is need to obtain a correct result (4325 -> 32633)
    """

    try:

        transformer_group = TransformerGroup(
            crs_from=f"epsg:{source_epsg_code}",
            crs_to=f"epsg:{dest_epsg_code}",
            area_of_interest=AreaOfInterest(*area_of_interest),
        )

        if not transformer_group.best_available:
            return False, "Best transformation is not available"

        #TODO: understand why you have to swap axis to obtain a correct result....
        proj_x_coords, proj_y_coords = transformer_group.transformers[0].transform(
            y_array,
            x_array
        )

        return True, (proj_x_coords, proj_y_coords)

    except Exception as e:

        return False, str(e)