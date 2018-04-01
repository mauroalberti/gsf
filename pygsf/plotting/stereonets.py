# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import mplstereonet as ms


from ..geotype_checks import *


def splot(data):
    """
    Plot geological data with matplotlib and mplstereonet
    """

    def params_gvect(gvects, upward):

        plunges = []
        bearings = []
        if upward:
            for gvect in gvects:
                plunge, bearing = gvect.mirror_horiz().pt
                plunges.append(plunge)
                bearings.append(bearing)
            symbol = "x"
        else:
            for gvect in gvects:
                plunge, bearing = gvect.pt
                plunges.append(plunge)
                bearings.append(bearing)
            symbol = "o"

        color = "blue"

        return np.array(plunges), np.array(bearings), symbol, color

    def params_gaxis(gaxes, upward):

        plunges = []
        bearings = []

        if upward:
            for gaxis in gaxes:
                plunge, bearing = gaxis.opposite().pt
                plunges.append(plunge)
                bearings.append(bearing)
        else:
            for gaxis in gaxes:
                plunge, bearing = gaxis.pt
                plunges.append(plunge)
                bearings.append(bearing)

        symbol = "s"
        color = "orange"

        return np.array(plunges), np.array(bearings), symbol, color

    if not isinstance(data, list):
        data = [data]

    gvects = list(filter(is_gvect, data))
    gaxes = list(filter(is_gaxis, data))

    gvects_notup = list(filter(is_not_upward, gvects))
    gvects_up = list(filter(is_upward, gvects))
    gaxes_notup = list(filter(is_not_upward, gaxes))
    gaxes_up = list(filter(is_upward, gaxes))
    # TODO: gplanes = filter(is_gplane, data)
    # TODO: slickens = filter(is_slickln, data)
    # TODO: faultslicks = filter(is_faultslck, data)
    # TODO: ptbaxes = filter(is_ptbaxes, data)

    if not any([gvects, gaxes]):
        return False

    gvects_up_params = params_gvect(gvects_up, upward=True) if gvects_up else None
    gvects_notup_params = params_gvect(gvects_notup, upward=False) if gvects_notup else None

    gaxes_up_params = params_gaxis(gaxes_up, upward=True) if gaxes_up else None
    gaxes_notup_params = params_gaxis(gaxes_notup, upward=False) if gaxes_notup else None

    datasets_1d = [
        gvects_up_params,
        gvects_notup_params,
        gaxes_up_params,
        gaxes_notup_params]

    fig, ax = ms.subplots()

    for data in datasets_1d:
        if data:
            plunges, bearings, symbol, color = data
            ax.line(
                plunges,
                bearings,
                marker=symbol,
                markerfacecolor=color,
                markeredgecolor=color)

    ax.grid()
    plt.show()


def splot_gvect(gvect):

    if gvect.is_upward:
        plunge, bearing = gvect.mirror_horiz().pt
        symbol = "x"
    else:
        plunge, bearing = gvect.pt
        symbol = "o"

    fig, ax = ms.subplots()

    ax.line(plunge, bearing, marker=symbol)
    ax.grid()
    plt.show()


def splot_gaxis(gaxis):

    pass


def splot_gplane(gplane):

    pass


def splot_slickln(slickln):
    pass


def splot_faultslick(faultslick):

    pass


class PlotException(Exception):

    pass

