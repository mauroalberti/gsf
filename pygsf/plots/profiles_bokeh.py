
from warnings import warn
from functools import singledispatch

from bokeh.plotting import figure, output_notebook, show
from bokeh.models import Range1d

from pygsf.profiles.geoprofiles import *


default_width = 18.5
default_height = 10.5


@singledispatch
def plot(
    obj,
    **kargs
) -> Optional[figure]:
    """

    :param obj:
    :param kargs:
    :return:
    """

    fig = kargs.get("fig", None)
    aspect = kargs.get("aspect", 1)
    width = kargs.get("width", default_width)
    height = kargs.get("height", default_height)

    if fig is None:

        output_notebook()
        fig = figure()

    show(fig)

    return fig


@plot.register(XYArrayPair)
def _(
    xyarrays: XYArrayPair,
    **kargs
) -> Optional[figure]:

    fig = kargs.get("fig", None)

    if fig is None:

        output_notebook()
        fig = figure()

    fig.match_aspect = True

    fig.line(
        xyarrays.x_arr(),
        xyarrays.y_arr(),
        line_width=0.75,
    )

    show(fig)

    return fig
