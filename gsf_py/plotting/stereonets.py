# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import mplstereonet as ms


from ..faults import *


def splot(data):
    """
    Plot data with matplotlib and mplstereonet
    """

    if isinstance(data, GVect):
        splot_gvect(data)
    elif isinstance(data, GAxis):
        splot_gaxis(data)
    elif isinstance(data, GPlane):
        splot_gplane(data)
    elif isinstance(data, Slickenline):
        splot_slickln(data)
    elif isinstance(data, FaultSlick):
        splot_faultslick(data)
    else:
        raise PlotException("Not available plot method for object type {}".format(type(data)))


def splot_gvect(data):

    fig, ax = ms.subplots()
    plunge, bearing = data.pt
    ax.line(plunge, bearing)
    ax.grid()
    plt.show()


def splot_gaxis(data):

    pass


def splot_gplane(data):

    pass


def splot_slickln(data):
    pass


def splot_faultslick(data):

    pass


class PlotException(Exception):

    pass

