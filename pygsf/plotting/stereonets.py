# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import mplstereonet as ms


from ..geotype_checks import *
from ..conversions.conversions import *


default_upward_marker = "x"
default_downward_marker = "o"

default_marker_color = "blue"


def splot(data, force=''):
    """
    Plot geological data with matplotlib and mplstereonet
    """

    def params_gvect(gvect, kwargs):

        if (force == 'lower' and gvect.is_upward) or \
           (force == 'upper' and gvect.is_downward):
            plot_gvect = gvect.opposite()
        else:
            plot_gvect = gvect

        if force == 'lower':
            default_marker = default_downward_marker
        elif force == 'upper':
            default_marker = default_upward_marker
        elif plot_gvect.is_upward:
            default_marker = default_upward_marker
        else:
            default_marker = default_downward_marker

        if plot_gvect.is_upward:  # apparently mplstereonet does not handle negative plunges
            plot_gvect = plot_gvect.mirror_horiz()

        plunge, bearing = plot_gvect.pt
        symbol = kwargs.get("m", default_marker)
        color = kwargs.get("c", default_marker_color)

        return plunge, bearing, symbol, color

    def params_gvects(gvects, upward):

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

    def params_gaxes(gaxes, upward):

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

    def params_gplanes(planes):

        strikes = []
        dips = []

        for gplane in gplanes:

            strike, dip = gplane.srda
            strikes.append(strike)
            dips.append(dip)

        color = "blue"

        return np.array(strikes), np.array(dips), color

    if force not in ('', 'upper', 'lower'):
        raise PlotException("Force parameter not valid")

    if not isinstance(data, list):
        data = [data]

    fig, ax = ms.subplots()

    for rec in data:

        if isinstance(rec, tuple):
            if isinstance(rec[-1], str):
                params = rec[-1]
                objs = rec[:-1]
            else:
                params = None
                objs = rec
        else:
            objs = [rec]
            params = None

        if params:
            kwargs = string2dict(params)
        else:
            kwargs = dict()

        for obj in objs:
            if is_gvect(obj):
                plunge, bearing, symbol, color = params_gvect(obj, kwargs)
                ax.line(
                    plunge,
                    bearing,
                    marker=symbol,
                    markerfacecolor=color,
                    markeredgecolor=color)


    """

    gvects = list(filter(is_gvect, data))
    gaxes = list(filter(is_gaxis, data))

    gvects_notup = list(filter(is_not_upward, gvects))
    gvects_up = list(filter(is_upward, gvects))
    gaxes_notup = list(filter(is_not_upward, gaxes))
    gaxes_up = list(filter(is_upward, gaxes))

    gplanes = filter(is_gplane, data)

    # TODO: slickens = filter(is_slickln, data)
    # TODO: faultslicks = filter(is_faultslck, data)
    # TODO: ptbaxes = filter(is_ptbaxes, data)

    if not any([gvects, gaxes, gplanes]):
        return False

    gvects_up_params = params_gvects(gvects_up, upward=True) if gvects_up else None
    gvects_notup_params = params_gvects(gvects_notup, upward=False) if gvects_notup else None

    gaxes_up_params = params_gaxes(gaxes_up, upward=True) if gaxes_up else None
    gaxes_notup_params = params_gaxes(gaxes_notup, upward=False) if gaxes_notup else None

    gplanes_params = params_gplanes(gplanes) if gplanes else None

    datasets_1d = [
        gvects_up_params,
        gvects_notup_params,
        gaxes_up_params,
        gaxes_notup_params]

    datasets_2d = [
        gplanes_params]
    """

    """
    for data in datasets_1d:
        if data:
            plunges, bearings, symbol, color = data
            ax.line(
                plunges,
                bearings,
                marker=symbol,
                markerfacecolor=color,
                markeredgecolor=color)
    """

    """
    for data in datasets_2d:
        if data:
            strikes, dips, color = data
            ax.plane(
                strikes,
                dips,
                color=color)
    """

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

