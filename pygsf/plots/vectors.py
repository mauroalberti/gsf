from matplotlib.figure import Figure

from pygsf.geometries.shapes.space3d import Line3D


def plot_line(
    fig: Figure,
    line: Line3D
) -> Figure:
    """
    Plot a line.

    :param fig: the figure in which to plot the line.
    :type fig: Figure.
    :param line: the line to plot.
    :type line: Line.
    :return: the input Figure instance.
    :rtype: Figure.
    """

    fig.gca().plot(line.x_list(), line.y_list(), '-')

    return fig