
min_epsg_crs_code = 2000  # checked 2019-06-14 in EPSG database


class Crs(object):
    """
    CRS class.
    Currently it is basic, just to manage simple comparisons and validity checks.

    """

    @staticmethod
    def undefined():

        return

    def __init__(self, epsg: int = -1):

        self._epsg = int(epsg)

    def epsg(self):

        return self._epsg

    def valid(self):

        return self.epsg() >= min_epsg_crs_code

    def __repr__(self):

        return "EPSG:{}".format(self.epsg())

    def __eq__(self, another):

        if not self.valid() or not another.valid():
            return False
        else:
            return self.epsg() == another.epsg()


