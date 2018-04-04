# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import mplstereonet as ms


from ..geotype_checks import *
from ..conversions.conversions import *


default_gvect_marker_upward_symbol = "x"
default_gvect_marker_downward_symbol = "o"
default_gvect_marker_color = "blue"

default_gaxis_marker_upward_symbol = "+"
default_gaxis_marker_downward_symbol = "s"
default_gaxis_marker_color = "orange"

default_gplane_downward_linestyle = "solid"
default_gplane_upward_linestyle = "dashed"
default_gplane_line_color = "blue"


def splot(data, force=''):
    """
    Plot geological data with matplotlib and mplstereonet
    """

    def params_gvect(gvect, kwargs, force_emisphere):

        if force_emisphere == 'lower':
            default_marker = default_gvect_marker_downward_symbol
            if gvect.is_upward:
                plot_gvect = gvect.opposite()
            else:
                plot_gvect = gvect
        elif force_emisphere == 'upper':
            default_marker = default_gvect_marker_upward_symbol
            if gvect.is_downward:
                plot_gvect = gvect.opposite()
            else:
                plot_gvect = gvect
        elif not force_emisphere:
            plot_gvect = gvect
            default_marker = default_gvect_marker_downward_symbol if not plot_gvect.is_upward else default_gvect_marker_upward_symbol
        else:
            raise PlotException("Invalid force emisphere parameter")

        if plot_gvect.is_upward:  # apparently mplstereonet does not handle negative plunges
            plot_gvect = plot_gvect.mirror_horiz()

        plunge, bearing = plot_gvect.pt
        symbol = kwargs.get("m", default_marker)
        color = kwargs.get("c", default_gvect_marker_color)

        return plunge, bearing, symbol, color

    def params_gaxis(gaxis, kwargs, force_emisphere):

        if (not force_emisphere) or (force_emisphere == 'lower'):
            default_marker = default_gaxis_marker_downward_symbol
            if gaxis.is_upward:
                plot_gaxis = gaxis.opposite()
            else:
                plot_gaxis = gaxis
        elif force_emisphere == 'upper':
            default_marker = default_gaxis_marker_upward_symbol
            if gaxis.is_downward:
                plot_gaxis = gaxis.opposite()
            else:
                plot_gaxis = gaxis
        else:
            raise PlotException("Invalid force emisphere parameter")

        if plot_gaxis.is_upward:  # apparently mplstereonet does not handle negative plunges
            plot_gaxis = plot_gaxis.mirror_horiz()

        plunge, bearing = plot_gaxis.pt
        symbol = kwargs.get("m", default_marker)
        color = kwargs.get("c", default_gaxis_marker_color)

        return plunge, bearing, symbol, color

    def params_gplane(gplane, kwargs, force_emisphere):

        if (not force_emisphere) or (force_emisphere == 'lower'):
            default_line_style = default_gplane_downward_linestyle
            plot_gplane = gplane
        elif force_emisphere == 'upper':
            default_line_style = default_gplane_upward_linestyle
            plot_gplane = gplane.mirror_vertical()
        else:
            raise PlotException("Invalid force emisphere parameter")

        strike, dip = plot_gplane.srda

        line_style = kwargs.get("m", default_line_style)
        color = kwargs.get("c", default_gplane_line_color)

        return strike, dip, line_style, color

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
                plunge, bearing, symbol, color = params_gvect(obj, kwargs, force_emisphere=force)
                ax.line(
                    plunge,
                    bearing,
                    marker=symbol,
                    markerfacecolor=color,
                    markeredgecolor=color)
            elif is_gaxis(obj):
                plunge, bearing, symbol, color = params_gaxis(obj, kwargs, force_emisphere=force)
                ax.line(
                    plunge,
                    bearing,
                    marker=symbol,
                    markerfacecolor=color,
                    markeredgecolor=color)
            elif is_gplane(obj):
                strike, dip, linestyle, color = params_gplane(obj, kwargs, force_emisphere=force)
                ax.plane(
                    strike,
                    dip,
                    linestyle=linestyle,
                    color=color)


    """

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

