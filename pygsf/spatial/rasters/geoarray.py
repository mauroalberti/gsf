# -*- coding: utf-8 -*-

from typing import Tuple, Optional, List

import numpy as np
from numpy import array

from ...defaults.types import Number

from ...mathematics.arrays import array_bilin_interp

from ..vectorial.geometries import Point
from ..exceptions import GeoArrayIOException

from ..projections.crs import Crs
from .geotransform import xyGeogrToijPix, gtToxyCellCenters, GeoTransform, ijPixToxyGeogr
from .fields import magnitude, orients_d, divergence, curl_module, magn_grads, magn_grad_along_flowlines


def ijPixToijArray(i_pix: Number, j_pix: Number) -> Tuple[Number, Number]:
    """
    Converts from pixel (geotransform-derived) to array indices.

    :param i_pix: the geotransform i value.
    :type i_pix: Number.
    :param j_pix: the geotransform j value.
    :type j_pix: Number.
    :return: the array-equivalent i and j indices.
    :rtype: a tuple of two numbers.

    Examples:
      >>> ijPixToijArray(0, 0)
      (-0.5, -0.5)
      >>> ijPixToijArray(0.5, 0.5)
      (0.0, 0.0)
      >>> ijPixToijArray(0.5, 1.5)
      (0.0, 1.0)
    """

    return i_pix - 0.5, j_pix - 0.5


def ijArrToijPix(i_arr: Number, j_arr: Number) -> Tuple[Number, Number]:
    """
    Converts from array indices to geotransform-related pixel indices.

    :param i_arr: the array i value.
    :type i_arr: Number.
    :param j_arr: the array j value.
    :type j_arr: Number.
    :return: the geotransform-equivalent i and j indices.
    :rtype: a tuple of two numbers.

    Examples:
      >>> ijArrToijPix(0, 0)
      (0.5, 0.5)
      >>> ijArrToijPix(0.5, 0.5)
      (1.0, 1.0)
      >>> ijArrToijPix(1.5, 0.5)
      (2.0, 1.0)
    """

    return i_arr + 0.5, j_arr + 0.5


class GeoArray(object):
    """
    GeoArray class.
    Stores and process georeferenced raster data.
    """

    def __init__(self, inGeotransform: GeoTransform, epsg_cd: int = -1, inLevels: Optional[List[np.ndarray]] = None) -> None:
        """
        GeoArray class constructor.

        :param  inGeotransform:  the geotransform
        :type  inGeotransform:  GeoTransform.
        :param epsg_cd: the projection EPSG code.
        :type epsg_cd: int
        :param  inLevels:  the nd-array storing the data.
        :type  inLevels:  np.array.

        :return:  None.

        Examples:
        """

        self._gt = inGeotransform
        self._crs = Crs(epsg_cd)
        if inLevels is None:
            self._levels = []
        else:
            self._levels = inLevels

    def geotransform(self):
        """
        Returns geotransform.

        :return: the geotransform.
        :rtype: GeoTransform.
        """

        return self._gt

    def crs(self):
        """
        Return the geoarray crs.

        :return: the crs.
        :rtype: Crs.
        """

        return self._crs

    def epsg(self) -> int:
        """
        Return the geoarray crs EPSG code.

        :return: the crs EPSG  code.
        :rtype: int.
        """

        return self._crs.epsg()

    def define_epsg(self, epsg_cd: int):
        """
        Overwrite the geoarray EPSG code.

        :return:
        """

        if not isinstance(epsg_cd, int):
            raise Exception("Provided EPSG code must be integer")

        self._crs = Crs(epsg_cd)

    def __repr__(self) -> str:
        """
        Represents a GeoArray instance as a shortened text.

        :return: a textual shortened representation of a GeoArray instance.
        :rtype: basestring.
        """

        num_bands = self.levels_num
        epsg_code = self.epsg()
        bands_txt = ""
        for band_ndx in range(num_bands):
            band = self.level(level_ndx=band_ndx)
            rows, cols = band.shape
            min, max = band.min(), band.max()
            bands_txt += "\nBand {}: {} rows x {} cols; min: {},  max: {}".format(band_ndx+1, rows, cols, min, max)

        txt = "GeoArray with {} band(s) - CRS: EPSG: {}\n{}".format(num_bands, epsg_code, bands_txt)

        return txt

    @property
    def src_cellsize_j(self) -> float:
        """
        Get the cell size of the geoarray in the x direction.

        :return: cell size in the x (j) direction.
        :rtype: float.

        Examples:
        """

        return abs(self._gt.pixWidth)

    @property
    def src_cellsize_i(self) -> float:
        """
        Get the cell size of the geoarray in the y direction.

        :return: cell size in the y (-i) direction.
        :rtype: float.

        Examples:
        """

        return abs(self._gt.pixHeight)

    @property
    def levels_num(self) -> int:
        """
        Returns the number of levels (dimensions) of the geoarray.

        :return: number of levels.
        :rtype: int.

        Examples:
          >>> gt = GeoTransform(0, 0, 10, 10)
          >>> GeoArray(gt, -1, [array([[1, 2], [3, 4]])]).levels_num
          1
          >>> GeoArray(gt, -1, [array([[1, 2], [3, 4]]), np.ones((4, 3, 2))]).levels_num
          2
        """

        return len(self._levels)

    def level(self, level_ndx: int=0):
        """
        Return the array corresponding to the requested level
        if existing else None.

        :param level_ndx: the index of the requested level.
        :type level_ndx: int.
        :return: the array or None.
        :rtype: optional array.

        Examples:
        """

        if 0 <= level_ndx < self.levels_num:
            return self._levels[level_ndx]
        else:
            return None

    def level_shape(self, level_ndx: int=0) -> Optional[Tuple[int, int]]:
        """
        Returns the shape (num. rows and num. columns) of the considered level grid.

        :param level_ndx: index of the level (grid) to consider.
        :type level_ndx: int.
        :return: number of rows and columns of the specific grid.
        :rtype: optional tuple of two int values.

        Examples:
          >>> gt = GeoTransform(0, 0, 10, 10)
          >>> GeoArray(gt, -1, [array([[1, 2], [3, 4]])]).level_shape()
          (2, 2)
          >>> GeoArray(gt, -1, [array([[1, 2], [3, 4]]), np.ones((4, 3, 2))]).level_shape(1)
          (4, 3, 2)
        """

        if 0 <= level_ndx < self.levels_num:
            return self._levels[level_ndx].shape
        else:
            return None

    def level_llc(self, level_ndx: int = 0) -> Optional[Tuple[int, int]]:
        """
        Deprecated. Use "band_corners_pixcoords" instead.

        Returns the coordinates of the lower-left corner.

        :param level_ndx: index of the level (grid) to consider.
        :type level_ndx: int.
        :return: x and y values of the lower-left corner of the specific grid.
        :rtype: optional tuple of two int values.

        Examples:
        """

        shape = self.level_shape(level_ndx)
        if not shape:
            return None

        llc_i_pix, llc_j_pix = shape[0], 0

        return self.ijPixToxy(llc_i_pix, llc_j_pix)

    def band_corners_pixcoords(self, level_ndx: int = 0) -> \
            Tuple[Tuple[float, float], Tuple[float, float], Tuple[float, float], Tuple[float, float]]:
        """
        Returns the pixel coordinates of the top-left, top-right, bottom-right and bottom-left band corners.

        :param level_ndx: index of the level (grid) to consider.
        :type level_ndx: int.
        :return: pixel coordinates of the top-left, top-right, bottom-right and bottom-left band corners.
        :rtype: four tuples of float pairs.

        Examples:
          >>> gt = GeoTransform(0, 0, 10, 10)
          >>> ga = GeoArray(gt, -1, [array([[1, 2, 3], [4, 5, 6]])])
          >>> ga.band_corners_pixcoords()
          ((0.0, 0.0), (0.0, 3.0), (2.0, 3.0), (2.0, 0.0))
        """

        shape = self.level_shape(level_ndx)
        num_rows, num_cols = shape

        top_left_ijpix = (0.0, 0.0)
        top_right_ijpix = (0.0, float(num_cols))
        btm_right_ijpix = (float(num_rows), float(num_cols))
        btm_left_ijpix = (float(num_rows), 0.0)

        return top_left_ijpix, top_right_ijpix, btm_right_ijpix, btm_left_ijpix

    def band_corners_geogcoords(self, level_ndx: int = 0) -> \
            Tuple[Tuple[float, float], Tuple[float, float], Tuple[float, float], Tuple[float, float]]:
        """
        Returns the geographic coordinates of the top-left, top-right, bottom-right and bottom-left band corners.

        :param level_ndx: index of the level (grid) to consider.
        :type level_ndx: int.
        :return: geographic coordinates of the top-left, top-right, bottom-right and bottom-left band corners.
        :rtype: four tuples of float pairs.

        Examples:
          >>> gt = GeoTransform(1500, 3000, 10, 10)
          >>> ga = GeoArray(gt, -1, [array([[1, 2, 3], [4, 5, 6]])])
          >>> ga.band_corners_geogcoords()
          ((1500.0, 3000.0), (1530.0, 3000.0), (1530.0, 2980.0), (1500.0, 2980.0))
        """

        top_left_ijpix, top_right_ijpix, btm_right_ijpix, btm_left_ijpix = self.band_corners_pixcoords(level_ndx=level_ndx)

        top_left_geogcoord = self.ijPixToxy(*top_left_ijpix)
        top_right_geogcoord = self.ijPixToxy(*top_right_ijpix)
        btm_right_geogcoord = self.ijPixToxy(*btm_right_ijpix)
        btm_left_geogcoord = self.ijPixToxy(*btm_left_ijpix)

        return top_left_geogcoord, top_right_geogcoord, btm_right_geogcoord, btm_left_geogcoord

    def xyToijArr(self, x: Number, y: Number) -> Tuple[Number, Number]:
        """
        Converts from geographic to array coordinates.

        :param x: x geographic component.
        :type x: Number.
        :param y: y geographic component.
        :type y: Number.
        :return: i and j values referred to array.
        :type: tuple of two float values.

        Examples:
        """

        return ijPixToijArray(*xyGeogrToijPix(self._gt, x, y))

    def xyToijPix(self, x: Number, y: Number) -> Tuple[Number, Number]:
        """
        Converts from geographic to pixel coordinates.

        :param x: x geographic component
        :type x: Number
        :param y: y geographic component
        :type y: Number
        :return: i and j values referred to grid.
        :type: tuple of two float values

        Examples:
        """

        return xyGeogrToijPix(self._gt, x, y)

    def ijArrToxy(self, i: Number, j: Number) -> Tuple[Number, Number]:
        """
        Converts from array indices to geographic coordinates.

        :param i: i array component.
        :type i: Number.
        :param j: j array component.
        :type j: Number.
        :return: x and y geographic coordinates.
        :type: tuple of two float values.

        Examples:
        """

        i_pix, j_pix = ijArrToijPix(i, j)

        return ijPixToxyGeogr(self._gt, i_pix, j_pix)

    def ijPixToxy(self, i: Number, j: Number) -> Tuple[Number, Number]:
        """
        Converts from grid indices to geographic coordinates.

        :param i: i pixel component.
        :type i: Number.
        :param j: j pixel component.
        :type j: Number.
        :return: x and y geographic coordinates.
        :type: tuple of two float values.

        Examples:
        """

        return ijPixToxyGeogr(self._gt, i, j)

    @property
    def has_rotation(self) -> bool:
        """
        Determines if a geoarray has axis rotations defined.

        :return: true if there are rotations, false otherwise.
        :rtype: bool.

        Examples:
        """

        return self._gt.has_rotation

    def geotransf_cell_sizes(self) -> Tuple[float, float]:
        """
        Calculates the geotransformed cell sizes.

        :return: a pair of float values, representing the cell sizes in the j and i directions.
        """

        factor = 100

        start_pt = Point(*self.ijArrToxy(0, 0))
        end_pt_j = Point(*self.ijArrToxy(0, factor))
        end_pt_i = Point(*self.ijArrToxy(factor, 0))

        return end_pt_j.dist2DWith(start_pt)/factor, end_pt_i.dist2DWith(start_pt)/factor

    def xy(self, level_ndx: int=0) -> Optional[Tuple[np.ndarray, np.ndarray]]:
        """
        Returns the two arrays storing respectively the x and the y coordinates
        of the grid cell centers for the chosen level (default is first level).

        :param level_ndx: the index of the
        :return: two arrays storing the geographical coordinates of the grid centers.
        :rtype: tuple made up by two float arrays.

        Examples:
        """

        res = self.level_shape(level_ndx)

        if not res:
            return None
        else:
            num_rows, num_cols = res
            return gtToxyCellCenters(self._gt, num_rows, num_cols)

    def interpolate_bilinear(self, x: Number, y: Number, level_ndx=0) -> Optional[float]:
        """
        Interpolate the z value at a point, given its geographic coordinates.
        Interpolation method: bilinear.

        :param x: x geographic coordinate.
        :type x: Number.
        :param y: y geographic coordinate.
        :type y: Number.
        :param level_ndx: the index of the used array.
        :type level_ndx: int.
        :return: the interpolated z value.
        :rtype: optional float.

        Examples:
        """

        i, j = self.xyToijArr(x, y)

        return array_bilin_interp(self._levels[level_ndx], i, j)

    def interpolate_bilinear_point(self, pt: Point, level_ndx=0) -> Optional[Point]:
        """
        Interpolate the z value at a point, returning a Point with elevation extracted from the DEM.
        Interpolation method: bilinear.

        :param pt: the positional point.
        :type pt: Point.
        :param level_ndx: the index of the used array.
        :type level_ndx: int.
        :return: a point with the same x-y position of the input point and with z equal to the interpolated z value.
        :rtype: optional Point.

        Examples:
        """

        x, y, epsg_cd = pt.x, pt.y, pt.epsg()

        z = self.interpolate_bilinear(x=x, y=y, level_ndx=level_ndx)

        if z:
            return Point(x, y, z, epsg_cd=epsg_cd)
        else:
            return None

    def magnitude_field(self, ndx_fx=0, ndx_fy=1) -> 'GeoArray':
        """
        Calculates magnitude field as a geoarray.

        :param ndx_fx: index of x field.
        :type ndx_fx: integer.
        :param ndx_fy: index of y field.
        :type ndx_fy: integer.
        :return: a geoarray storing the magnitude field.
        :rtype: GeoArray.

        Examples:
        """

        magn = magnitude(
            fld_x=self._levels[ndx_fx],
            fld_y=self._levels[ndx_fy])

        return GeoArray(
            inGeotransform=self._gt,
            epsg_cd=self.epsg(),
            inLevels=[magn]
        )

    def orientations(self, ndx_fx=0, ndx_fy=1) -> 'GeoArray':
        """
        Calculates orientations field as a geoarray.

        :param ndx_fx: index of x field.
        :type ndx_fx: integer.
        :param ndx_fy: index of y field.
        :type ndx_fy: integer.
        :return: a geoarray storing the orientation field.
        :rtype: GeoArray.

        Examples:
        """

        orient = orients_d(
            fld_x=self._levels[ndx_fx],
            fld_y=self._levels[ndx_fy])

        return GeoArray(
            inGeotransform=self._gt,
            epsg_cd=self._crs,
            inLevels=[orient]
        )

    def divergence_2D(self, ndx_fx=0, ndx_fy=1) -> 'GeoArray':
        """
        Calculates divergence of a 2D field as a geoarray.

        :param ndx_fx: index of x field.
        :type ndx_fx: integer.
        :param ndx_fy: index of y field.
        :type ndx_fy: integer.
        :return: a geoarray storing the divergence field.
        :rtype: GeoArray.

        Examples:
        """

        div = divergence(
            fld_x=self._levels[ndx_fx],
            fld_y=self._levels[ndx_fy],
            cell_size_x=self.src_cellsize_j,
            cell_size_y=self.src_cellsize_i)

        return GeoArray(
            inGeotransform=self._gt,
            epsg_cd=self._crs,
            inLevels=[div]
        )

    def curl_module(self, ndx_fx=0, ndx_fy=1) -> 'GeoArray':
        """
        Calculates curl module of a 2D field as a geoarray.

        :param ndx_fx: index of x field.
        :type ndx_fx: integer.
        :param ndx_fy: index of y field.
        :type ndx_fy: integer.
        :return: a geoarray storing the curl module field.
        :rtype: GeoArray.

        Examples:
        """

        curl_m = curl_module(
            fld_x=self._levels[ndx_fx],
            fld_y=self._levels[ndx_fy],
            cell_size_x=self.src_cellsize_j,
            cell_size_y=self.src_cellsize_i)

        return GeoArray(
            inGeotransform=self._gt,
            epsg_cd=self._crs,
            inLevels=[curl_m])

    def magnitude_grads(self, axis: str= '', ndx_fx: int=0, ndx_fy: int=1) -> 'GeoArray':
        """
        Calculates the magnitude gradient along the x, y axis or both, of a 2D field as a geoarray.

        :param axis: axis along wich to calculate the gradient, 'x' or 'y', or '' (predefined) for both x and y.
        :type axis: str.
        :param ndx_fx: index of x field.
        :type ndx_fx: integer.
        :param ndx_fy: index of y field.
        :type ndx_fy: integer.
        :return: a geoarray storing the magnitude gradient along the x, y axis (or both) field.
        :rtype: GeoArray.
        :raises: GeoArrayIOException.

        Examples:
        """

        if axis == 'x':
            cell_sizes = [self.src_cellsize_j]
        elif axis == 'y':
            cell_sizes = [self.src_cellsize_i]
        elif axis == '':
            cell_sizes = [self.src_cellsize_j, self.src_cellsize_i]
        else:
            raise GeoArrayIOException("Axis must be 'x' or 'y. '{}' given".format(axis))

        magnitude_gradients = magn_grads(
            fld_x=self._levels[ndx_fx],
            fld_y=self._levels[ndx_fy],
            dir_cell_sizes=cell_sizes,
            axis=axis)

        return GeoArray(
            inGeotransform=self._gt,
            epsg_cd=self._crs,
            inLevels=magnitude_gradients)

    def grad_flowlines(self, ndx_fx: int=0, ndx_fy: int=1) -> 'GeoArray':
        """
        Calculates gradient along flow lines.

        :param ndx_fx: index of x field.
        :type ndx_fx: integer.
        :param ndx_fy: index of y field.
        :type ndx_fy: integer.
        :return: a geoarray storing the flowline gradient field
        :rtype: GeoArray
        """

        flowln_grad = magn_grad_along_flowlines(
            fld_x=self._levels[ndx_fx],
            fld_y=self._levels[ndx_fy],
            cell_size_x=self.src_cellsize_j,
            cell_size_y=self.src_cellsize_i)

        return GeoArray(
            inGeotransform=self._gt,
            epsg_cd=self._crs,
            inLevels=[flowln_grad])


if __name__ == "__main__":

    import doctest
    doctest.testmod()

