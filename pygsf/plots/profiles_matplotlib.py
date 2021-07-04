
from warnings import warn
from functools import singledispatch

from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from matplotlib import gridspec

from pygsf.profiles.geoprofiles import *


colors_addit = [
    "blue",
    "darkseagreen",
    "darkgoldenrod",
    "darkviolet",
    "hotpink",
    "powderblue",
    "yellowgreen",
    "palevioletred",
    "seagreen",
    "darkturquoise",
    "beige",
    "darkkhaki",
    "red",
    "yellow",
    "magenta",
    "cyan",
    "chartreuse"
]

z_padding = 0.2

default_width = 18.5
default_height = 10.5


@singledispatch
def plot(
    obj,
    **kargs
) -> Optional[Figure]:
    """

    :param obj:
    :param kargs:
    :return:
    """

    fig = kargs.pop("fig", None)
    aspect = kargs.pop("aspect", 1)
    width = kargs.pop("width", default_width)
    height = kargs.pop("height", default_height)

    if fig is None:

        fig, ax = plt.subplots()
        fig.set_size_inches(width, height)

        ax.set_aspect(aspect)

    else:

        plt.gca()

    return fig


@plot.register(XYArrayPair)
def _(
    xyarrays: XYArrayPair,
    **kargs
) -> Optional[Figure]:

    fig = kargs.pop("fig", None)
    z_min = kargs.pop("z_min", None)
    z_max = kargs.pop("z_max", None)
    aspect = kargs.pop("aspect", 1)
    width = kargs.pop("width", default_width)
    height = kargs.pop("height", default_height)
    grid = kargs.pop("grid", True)
    breaklines = kargs.pop("breaklines", True)
    breaklines_color = kargs.pop("breaklines_color", 'tan')
    breaklines_width = kargs.pop("breaklines_width", 0.9)
    breaklines_style = kargs.pop("breaklines_style", 'dotted')

    if z_min is None or z_max is None:
        z_range = xyarrays.y_max() - xyarrays.y_min()
        z_min = xyarrays.y_min() - z_padding * z_range
        z_max = xyarrays.y_max() + z_padding * z_range

    if np.isnan(z_min) or np.isnan(z_max):
        return

    if fig is None:

        fig = plt.figure()
        fig.set_size_inches(width, height)

    ax = fig.add_subplot()

    ax.set_aspect(aspect)

    if z_min is not None or z_max is not None:
        ax.set_ylim([z_min, z_max])

    if grid:
        ax.grid(True, linestyle='-', color='0.90')

    ax.plot(
        xyarrays.x_arr(),
        xyarrays.y_arr(),
        **kargs
    )

    if breaklines:
        bottom, top = ax.get_ylim()
        ax.vlines(
            xyarrays.x_breaks(),
            bottom,
            top,
            color=breaklines_color,
            linewidth=breaklines_width,
            linestyles=breaklines_style
        )

    return fig


@plot.register(GeoProfile)
def _(
    geoprofile: GeoProfile,
    **kargs
) -> Optional[Figure]:
    """
    Plot a single geological profile.

    :param geoprofile: the geoprofile to plot
    :type geoprofile: qygsf.geometries.shapes.collections.NamedLines
    :return: the figure.
    :rtype: Figure
    """

    fig = kargs.pop("fig", None)
    z_min = kargs.pop("z_min", None)
    z_max = kargs.pop("z_max", None)
    ndx = kargs.pop("ndx", 0)
    aspect = kargs.pop("aspect", 1)
    width = kargs.pop("width", default_width)
    height = kargs.pop("height", default_height)
    superposed = kargs.pop("superposed", False)
    num_subplots = kargs.pop("num_subplots", 1)
    spec = kargs.pop("spec", None)
    labels_add_orientdip = kargs.pop("labels_add_orientdip", None)
    labels_add_id = kargs.pop("labels_add_id", None)
    inters_color = kargs.pop("inters_color", None)
    inters_label = kargs.pop("inters_label", None)
    attitude_color = kargs.pop("attitude_color", "red")

    if z_min is None or z_max is None:
        z_range = geoprofile.z_max() - geoprofile.z_min()
        z_min = geoprofile.z_min() - z_padding * z_range
        z_max = geoprofile.z_max() + z_padding * z_range

    if np.isnan(z_min) or np.isnan(z_max):
        return

    if fig is None:

        fig = plt.figure()
        fig.set_size_inches(width, height)

    if superposed:
        ax = fig.add_axes(
            [0.1, 0.1, 0.8, 0.8]
        )
    elif spec is not None:
        ax = fig.add_subplot(
            spec[ndx, 0]
        )
    else:
        ax = fig.add_subplot()

    ax.set_aspect(aspect)

    if z_min is not None or z_max is not None:
        ax.set_ylim([z_min, z_max])

    if geoprofile.topo_profile:

        kargs['color'] = colors_addit[ndx % len(colors_addit)] if superposed else colors_addit[-1]

        ax.plot(
            geoprofile.topo_profile.x_arr(),
            geoprofile.topo_profile.y_arr(),
            **kargs
        )

    if geoprofile.profile_attitudes:

        attits = geoprofile.profile_attitudes

        section_length = geoprofile.length_2d()

        projected_z = [structural_attitude.z for structural_attitude in attits if
                       0.0 <= structural_attitude.s <= section_length]

        projected_s = [structural_attitude.s for structural_attitude in attits if
                       0.0 <= structural_attitude.s <= section_length]

        projected_ids = [structural_attitude.id for structural_attitude in attits if
                         0.0 <= structural_attitude.s <= section_length]

        axes = fig.gca()
        vertical_exaggeration = axes.get_aspect()

        axes.plot(
            projected_s,
            projected_z,
            'o',
            color=attitude_color
        )

        # plot segments representing structural data

        for structural_attitude in attits:
            if 0.0 <= structural_attitude.s <= section_length:

                structural_segment_s, structural_segment_z = structural_attitude.create_segment_for_plot(
                    section_length,
                    vertical_exaggeration)

                fig.gca().plot(
                    structural_segment_s,
                    structural_segment_z,
                    '-',
                    color=attitude_color
                )

        if labels_add_orientdip or labels_add_id:

            src_dip_dirs = [structural_attitude.src_dip_dir for structural_attitude in
                            attits if 0.0 <= structural_attitude.s <= section_length]
            src_dip_angs = [structural_attitude.src_dip_ang for structural_attitude in
                            attits if 0.0 <= structural_attitude.s <= section_length]

            for rec_id, src_dip_dir, src_dip_ang, s, z in zip(
                    projected_ids,
                    src_dip_dirs,
                    src_dip_angs,
                    projected_s,
                    projected_z):

                if labels_add_orientdip and labels_add_id:
                    label = "%s-%03d/%02d" % (rec_id, src_dip_dir, src_dip_ang)
                elif labels_add_id:
                    label = "%s" % rec_id
                elif labels_add_orientdip:
                    label = "%03d/%02d" % (src_dip_dir, src_dip_ang)
                else:
                    raise Exception(f"Unhandled case with {labels_add_orientdip} and {labels_add_id}")

                axes.annotate(label, (s + 15, z + 15))

    if geoprofile.lines_intersections:

        if not geoprofile.topo_profile:

            warn('Topographic profile is not defined, so intersections cannot be plotted')

        else:

            for ndx, intersection_element in enumerate(geoprofile.lines_intersections):

                intersection_id = intersection_element.id
                intersection_subparts = intersection_element.arrays

                for s_range in intersection_subparts:

                    s_start = s_range[0]
                    s_end = s_range[1] if len(s_range) > 1 else None
                    plot_symbol = '-o' if len(s_range) > 1 else 'o'

                    s_vals = geoprofile.topo_profile.x_subset(
                        s_start,
                        s_end
                    )

                    z_vals = geoprofile.topo_profile.ys_from_x_range(
                        s_start,
                        s_end
                    )

                    fig.gca().plot(
                        s_vals,
                        z_vals,
                        plot_symbol,
                        color=colors_addit[ndx]
                    )

                    if inters_label:

                        fig.gca().annotate(
                            f"{intersection_id}",
                            (s_vals[-1] + 20, z_vals[-1] + 25))

    if geoprofile.polygons_intersections:

        if not geoprofile.topo_profile:

            warn('Topographic profile is not defined, so intersections cannot be plotted')

        else:

            for ndx, intersection_element in enumerate(geoprofile.polygons_intersections):

                intersection_id = intersection_element.id
                intersection_subparts = intersection_element.arrays

                for s_range in intersection_subparts:

                    s_start = s_range[0]
                    s_end = s_range[1] if len(s_range) > 1 else None
                    plot_symbol = '-' if len(s_range) > 1 else 'o'

                    s_vals = geoprofile.topo_profile.x_subset(
                        s_start,
                        s_end
                    )

                    z_vals = geoprofile.topo_profile.ys_from_x_range(
                        s_start,
                        s_end
                    )

                    fig.gca().plot(
                        s_vals,
                        z_vals,
                        plot_symbol,
                        color=colors_addit[ndx]
                    )

                    if inters_label:

                        fig.gca().annotate(
                            f"{intersection_id}",
                            (s_vals[-1] + 20, z_vals[-1] + 25))

    return fig


@plot.register(GeoProfileSet)
def _(
    geoprofiles: GeoProfileSet,
    **kargs
) -> List[Optional[Figure]]:
    """
    Plot a set of geological profiles.

    :param geoprofiles: the geoprofiles to plot
    :type geoprofiles: GeoProfiles
    :return: the figures.
    :rtype: List[Figure]
    """

    fig = kargs.get("fig", None)
    plot_z_min = kargs.get("plot_z_min", None)
    plot_z_max = kargs.get("plot_z_max", None)
    ndx = kargs.get("ndx", 0)
    aspect = kargs.get("aspect", 1)
    width = kargs.get("width", default_width)
    height = kargs.get("height", default_height)
    superposed = kargs.get("superposed", False)
    num_subplots = kargs.get("num_subplots", 1)
    spec = kargs.get("spec", None)
    labels_add_orientdip = kargs.get("labels_add_orientdip", None)
    labels_add_id = kargs.get("labels_add_id", None)

    if plot_z_min is None or plot_z_max is None:

        z_range = geoprofiles.z_max() - geoprofiles.z_min()
        plot_z_min = geoprofiles.z_min() - z_padding * z_range
        plot_z_max = geoprofiles.z_max() + z_padding * z_range

    if np.isnan(plot_z_min) or np.isnan(plot_z_max):
        return []

    num_profiles = geoprofiles.num_profiles()

    if not superposed:
        fig = plt.figure(constrained_layout=True)
        spec = gridspec.GridSpec(
            ncols=1,
            nrows=num_profiles,
            figure=fig)
    else:
        fig = plt.figure()
        spec = None

    fig.set_size_inches(width, height)

    for ndx in range(num_profiles):

        geoprofile = geoprofiles.extract_geoprofile(ndx)

        fig = plot(
            geoprofile,
            plot_z_min=plot_z_min,
            plot_z_max=plot_z_max,
            fig=fig,
            ndx=ndx,
            superposed=superposed,
            num_subplots=num_profiles,
            spec=spec,
            labels_add_orientdip=labels_add_orientdip,
            labels_add_id=labels_add_id
        )

    return fig

