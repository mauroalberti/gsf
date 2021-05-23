
import os

from pygsf.geometries.shapes.space2d import Point2D

try:
    from osgeo import gdal
except ImportError:
    import gdal

from ...georeferenced.rasters import *
from ...geometries.grids.geotransform import *


GRID_NULL_VALUE = -99999  # should already be imported but it isn't


def read_raster(
        file_ref: Any
) -> Tuple[gdal.Dataset, Optional[GeoTransform], int, str]:
    """
    Read a raster layer.

    :param file_ref: the reference to the raster
    :type file_ref: Any
    :return: the dataset, its geotransform, the number of bands, the projection.
    :rtype: tuple made up by a gdal.Dataset instance, an optional Geotransform object, and int and a string.
    :raises: Exception

    Examples:
    """

    # open raster file and check operation success

    dataset = gdal.Open(file_ref, gdal.GA_ReadOnly)
    if not dataset:
        raise Exception("No input data open")

    # get raster descriptive infos

    gt = dataset.GetGeoTransform()
    if gt:
        geotransform = GeoTransform.fromGdalGt(gt)
    else:
        geotransform = None

    num_bands = dataset.RasterCount

    projection = dataset.GetProjection()

    return dataset, geotransform, num_bands, projection


def read_band(
        dataset: gdal.Dataset,
        bnd_ndx: int = 1
) -> Tuple[dict, 'np.array']:
    """
    Read data and metadata of a grids band based on GDAL.

    :param dataset: the source raster dataset
    :type dataset: gdal.Dataset
    :param bnd_ndx: the index of the band (starts from 1)
    :type bnd_ndx: int
    :return: the band parameters and the data values
    :rtype: dict of data parameters and values as a numpy.array
    :raises: Exception

    Examples:

    """

    band = dataset.GetRasterBand(bnd_ndx)
    data_type = gdal.GetDataTypeName(band.DataType)

    unit_type = band.GetUnitType()

    stats = band.GetStatistics(False, False)
    if stats is None:
        dStats = dict(
            min=None,
            max=None,
            mean=None,
            std_dev=None)
    else:
        dStats = dict(
            min=stats[0],
            max=stats[1],
            mean=stats[2],
            std_dev=stats[3])

    noDataVal = band.GetNoDataValue()

    nOverviews = band.GetOverviewCount()

    colorTable = band.GetRasterColorTable()

    if colorTable:
        nColTableEntries = colorTable.GetCount()
    else:
        nColTableEntries = 0

    # read data from band

    grid_values = band.ReadAsArray()
    if grid_values is None:
        raise Exception("Unable to read data from grids")

    # transform data into numpy array

    data = np.asarray(grid_values)

    # if nodatavalue exists, set null values to NaN in numpy array
    if noDataVal is not None and np.isfinite(noDataVal):
        data = np.where(abs(data - noDataVal) > 1e-10, data, np.NaN)

    band_params = dict(
        dataType=data_type,
        unitType=unit_type,
        stats=dStats,
        noData=noDataVal,
        numOverviews=nOverviews,
        numColorTableEntries=nColTableEntries)

    return band_params, data


def try_read_raster_band(
    raster_source: str,
    bnd_ndx: int = 1
) -> Tuple[bool, Union[str, Tuple[GeoTransform, str, Dict, 'np.array']]]:

    # get raster parameters and data
    try:
        dataset, geotransform, num_bands, projection = read_raster(raster_source)
    except (IOError, TypeError, Exception) as err:
        return False, "Exception with reading {}: {}".format(raster_source, err)

    band_params, data = read_band(dataset, bnd_ndx)

    return True, (geotransform, projection, band_params, data)


def try_write_esrigrid(
    geoarray: GeoArray,
    outgrid_fn: str,
    esri_nullvalue: numbers.Integral = GRID_NULL_VALUE,
    level_ndx: int = 0
) -> Tuple[bool, str]:
    """
    Writes ESRI ascii grid.
    
    :param geoarray: 
    :param outgrid_fn: 
    :param esri_nullvalue: 
    :param level_ndx: index of the level array to write.
    :type level_ndx: int.
    :return: success and descriptive message
    :rtype: tuple made up by a boolean and a string
    """
    
    outgrid_fn = str(outgrid_fn)

    # checking existence of output slope grid

    if os.path.exists(outgrid_fn):
        return False, "Output grid '{}' already exists".format(outgrid_fn)

    try:
        outputgrid = open(outgrid_fn, 'w')  # create the output ascii file
    except Exception:
        return False, "Unable to create output grid '{}'".format(outgrid_fn)

    if outputgrid is None:
        return False, "Unable to create output grid '{}'".format(outgrid_fn)

    if geoarray.has_rotation:
        return False, "Grid has axes rotations defined"

    cell_size_x = geoarray.src_cellsize_j
    cell_size_y = geoarray.src_cellsize_i

    if not areClose(cell_size_x, cell_size_y):
        return False, "Cell sizes in the x- and y- directions are not similar"

    arr = geoarray.level(level_ndx)
    if arr is None:
        return False, f"Array with index {level_ndx} does not exist"

    num_rows, num_cols = arr.shape
    llc_x, llc_y = geoarray.level_llc(level_ndx)

    # writes header of grid ascii file

    outputgrid.write("NCOLS %d\n" % num_cols)
    outputgrid.write("NROWS %d\n" % num_rows)
    outputgrid.write("XLLCORNER %.8f\n" % llc_x)
    outputgrid.write("YLLCORNER %.8f\n" % llc_y)
    outputgrid.write("CELLSIZE %.8f\n" % cell_size_x)
    outputgrid.write("NODATA_VALUE %f\n" % esri_nullvalue)

    esrigrid_outvalues = np.where(np.isnan(arr), esri_nullvalue, arr)

    # output of results

    for i in range(0, num_rows):
        for j in range(0, num_cols):
            outputgrid.write("%.8f " % (esrigrid_outvalues[i, j]))
        outputgrid.write("\n")

    outputgrid.close()

    return True, "Data saved in {}".format(outgrid_fn)


if __name__ == "__main__":

    import doctest
    doctest.testmod()


class GDALParameters(object):
    """
    Manage GDAL parameters from grids.

    """

    # class constructor
    def __init__(self):
        """
        Class constructor.

        @return:  generic-case GDAL parameters.
        """
        self._nodatavalue = None
        self._topleftX = None
        self._topleftY = None
        self._pixsizeEW = None
        self._pixsizeNS = None
        self._rows = None
        self._cols = None
        self._rotation_GT_2 = 0.0
        self._rotation_GT_4 = 0.0

    def s_noDataValue(self, nodataval):
        """
        Set raster no data value.

        @param  nodataval:  the raster no-data value.
        @type  nodataval:  None, otherwise number or string convertible to float.

        @return:  self.
        """

        try:
            self._nodatavalue = float(nodataval)
        except:
            self._nodatavalue = None

    def g_noDataValue(self):
        """
        Get raster no-data value.

        @return:  no-data value - float.
        """
        return self._nodatavalue

    # set property for no-data value
    noDataValue = property(g_noDataValue, s_noDataValue)

    def s_topLeftX(self, topleftX):
        """
        Set top-left corner x value of the raster.

        @param  topleftX:  the top-left corner x value, according to GDAL convention.
        @type  topleftX:  number or string convertible to float.

        @return:  self.
        """
        self._topleftX = float(topleftX)

    def g_topLeftX(self):
        """
        Get top-left corner x value of the raster.

        @return:  the top-left corner x value, according to GDAL convention - float.
        """
        return self._topleftX

    # set property for topleftX
    topLeftX = property(g_topLeftX, s_topLeftX)

    def s_topLeftY(self, topleftY):
        """
        Set top-left corner y value of the raster.

        @param  topleftY:  the top-left corner y value, according to GDAL convention.
        @type  topleftY:  number or string convertible to float.

        @return:  self.
        """
        self._topleftY = float(topleftY)

    def g_topLeftY(self):
        """
        Get top-left corner y value of the raster.

        @return:  the top-left corner y value, according to GDAL convention - float.
        """
        return self._topleftY

    # set property for topleftY
    topLeftY = property(g_topLeftY, s_topLeftY)

    def s_pixSizeEW(self, pixsizeEW):
        """
        Set East-West size of the raster cell.

        @param  pixsizeEW:  the top-left y value, according to GDAL convention.
        @type  pixsizeEW:  number or string convertible to float.

        @return:  self.
        """
        self._pixsizeEW = float(pixsizeEW)

    def g_pixSizeEW(self):
        """
        Get East-West size of the raster cell.

        @return:  the East-West size of the raster cell - float.
        """
        return self._pixsizeEW

    # set property for topleftY
    pixSizeEW = property(g_pixSizeEW, s_pixSizeEW)

    # pixsizeNS

    def s_pixSizeNS(self, pixsizeNS):
        """
        Set North-South size of the raster cell.

        @param  pixsizeNS:  the North-South size of the raster cell.
        @type  pixsizeNS:  number or string convertible to float.

        @return:  self.
        """
        self._pixsizeNS = float(pixsizeNS)

    def g_pixSizeNS(self):
        """
        Get North-South size of the raster cell.

        @return:  the North-South size of the raster cell - float.
        """
        return self._pixsizeNS

    # set property for topleftY
    pixSizeNS = property(g_pixSizeNS, s_pixSizeNS)

    def s_rows(self, rows):
        """
        Set row number.

        @param  rows:  the raster row number.
        @type  rows:  number or string convertible to int.

        @return:  self.
        """
        self._rows = int(rows)

    def g_rows(self):
        """
        Get row number.

        @return:  the raster row number - int.
        """
        return self._rows

    # set property for rows
    rows = property(g_rows, s_rows)

    def s_cols(self, cols):
        """
        Set column number.

        @param  cols:  the raster column number.
        @type  cols:  number or string convertible to int.

        @return:  self.
        """
        self._cols = int(cols)

    def g_cols(self):
        """
        Get column number.

        @return:  the raster column number - int.
        """
        return self._cols

    # set property for cols
    cols = property(g_cols, s_cols)

    def s_rotation_GT_2(self, rotation_GT_2):
        """
        Set rotation GT(2) (see GDAL documentation).

        @param  rotation_GT_2:  the raster rotation value GT(2).
        @type  rotation_GT_2:  number or string convertible to float.

        @return:  self.
        """
        self._rotation_GT_2 = float(rotation_GT_2)

    def g_rotation_GT_2(self):
        """
        Get rotation GT(2) (see GDAL documentation).

        @return:  the raster rotation value GT(2). - float.
        """
        return self._rotation_GT_2

    # set property for rotation_GT_2
    rotGT2 = property(g_rotation_GT_2, s_rotation_GT_2)

    def s_rotation_GT_4(self, rotation_GT_4):
        """
        Set rotation GT(4) (see GDAL documentation)

        @param  rotation_GT_4:  the raster rotation value GT(4).
        @type  rotation_GT_4:  number or string convertible to float.

        @return:  self.
        """
        self._rotation_GT_4 = float(rotation_GT_4)

    def g_rotation_GT_4(self):
        """
        Get rotation GT(4) (see GDAL documentation).

        @return:  the raster rotation value GT(4) - float.
        """
        return self._rotation_GT_4

    # set property for rotation_GT_4
    rotGT4 = property(g_rotation_GT_4, s_rotation_GT_4)

    def check_params(self, tolerance=1e-06):
        """
        Check absence of axis rotations or pixel size differences in the raster band.

        @param  tolerance:  the maximum threshold for both pixel N-S and E-W difference, or axis rotations.
        @type  tolerance:  float.

        @return:  None when successful, RasterParametersException when pixel differences or axis rotations.

        @raise: RasterParametersException - raster geometry incompatible with this module (i.e. different cell sizes or axis rotations).
        """
        # check if pixel size can be considered the same in the two axis directions
        if abs(abs(self._pixsizeEW) - abs(self._pixsizeNS)) / abs(self._pixsizeNS) > tolerance:
            raise Exception('Pixel sizes in x and y directions are different in raster')

            # check for the absence of axis rotations
        if abs(self._rotation_GT_2) > tolerance or abs(self._rotation_GT_4) > tolerance:
            raise Exception('There should be no axis rotation in raster')

        return

    def llcorner(self):
        """
        Creates a point at the lower-left corner of the raster.

        @return:  new Point instance.
        """
        return Point2D(self.topLeftX, self.topLeftY - abs(self.pixSizeNS) * self.rows)

    def trcorner(self):
        """
        Create a point at the top-right corner of the raster.

        @return:  new Point instance.
        """
        return Point2D(self.topLeftX + abs(self.pixSizeEW) * self.cols, self.topLeftY)

    def geo_equiv(self, other, tolerance=1.0e-6):
        """
        Checks if two grids are geographically equivalent.

        @param  other:  a grid to be compared with self.
        @type  other:  Grid instance.
        @param  tolerance:  the maximum threshold for pixel sizes, topLeftX or topLeftY differences.
        @type  tolerance:  float.

        @return:  Boolean.
        """
        if 2 * (self.topLeftX - other.topLeftX) / (self.topLeftX + other.topLeftX) > tolerance or \
                                        2 * (self.topLeftY - other.topLeftY) / (
                                    self.topLeftY + other.topLeftY) > tolerance or \
                                        2 * (abs(self.pixSizeEW) - abs(other.pixSizeEW)) / (
                                    abs(self.pixSizeEW) + abs(other.pixSizeEW)) > tolerance or \
                                        2 * (abs(self.pixSizeNS) - abs(other.pixSizeNS)) / (
                                    abs(self.pixSizeNS) + abs(other.pixSizeNS)) > tolerance or \
                        self.rows != other.rows or self.cols != other.cols:
            return False
        else:
            return True