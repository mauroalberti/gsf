# -*- coding: utf-8 -*-


class CRSCodeException(Exception):
    """
    Class for exceptions related to CRS codes.
    """

    pass


class RasterIOExceptions(Exception):
    """
    Class for rasters IO exception
    """

    pass


class RasterParametersExceptions(Exception):
    """
    Class for rasters IO exception
    """

    pass


class AnaliticSurfaceIOException(Exception):
    """
    Exception for Analytical Surfaces IO
    """

    pass


class AnaliticSurfaceCalcException(Exception):
    """
    Exception for Analytical Surface calculation.
    """

    pass


class VectorIOException(Exception):

    pass


class GeoArrayIOException(Exception):
    """
    Class for geoarray IO exception
    """

    pass
