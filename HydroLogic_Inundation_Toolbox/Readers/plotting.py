from typing import List

import contextily as ctx
import matplotlib.pyplot as plt
import rasterio
from matplotlib.axes import Axes
from matplotlib_scalebar.scalebar import ScaleBar
from rasterio.plot import show


def default_plot(dfs: List, variable: str, labels: List = None) -> None:
    """Basic plot function that takes a number of DataFrames and plots the 'variable' columns in one figure

    Args:
        dfs (List): list of dataframes (can be 1) that are to be plotted. Should all contain a column with variable
        variable (str): name of variable and column to be selected from dataframes
        labels (List): labels of plotted data frames, e.g. ["Measured", "Simulated"]

    Returns:
        None
    """
    plt.figure(figsize=(12, 4))
    for ix, df in enumerate(dfs):
        if labels is None:
            plt.plot(df[variable].dropna())
        else:
            plt.plot(df[variable].dropna(), label=labels[ix])

    if labels is not None:
        plt.legend()

    plt.gca().update(dict(title=r"Plot of: " + variable, xlabel="date", ylabel=variable))
    plt.grid()
    plt.show()


def raster_plot_with_context(
    raster_path: str,
    epsg: int = 28992,
    basemap: bool = True,
    clabel: str = None,
    cmap: str = None,
    title: str = None,
    vmin: float = None,
    vmax: float = None,
) -> plt.figure:
    """
    Function to create a plot from a tiff file

    Args:
        raster_path (str): location of tiff file
        epsg (int): EPSG identifier of coordinate reference system. Default is EPSG:28992, Dutch RDS
        basemap (bool): Whether to load a basemap or not
        clabel (str): label for colorbar
        cmap (str): string identifier for matplotlib colormap
        title (str): title of plot
        vmin (float): lower bound of colorbar
        vmax (float): upper bound of colorbar

    Returns:
        fig (plt.Figure): matplotlib figure object
        ax (plt.axis): matplotlib axis

    """

    kwargs = {}
    if cmap is not None:
        kwargs["cmap"] = cmap
    if title is not None:
        kwargs["title"] = title
    if vmin is not None:
        kwargs["vmin"] = vmin
    if vmax is not None:
        kwargs["vmax"] = vmax

    figsize = [12, 6]
    dpi = 150

    fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=dpi)

    src = rasterio.open(raster_path)

    show(source=src, ax=ax, zorder=1, **kwargs)
    ax.invert_yaxis()
    cbar = plt.colorbar(ax.get_images()[0], ax=ax)

    if clabel is not None:
        cbar.set_label(clabel)

    ax.set(xlabel="x", ylabel="y")

    if basemap:
        ctx.add_basemap(ax, crs=epsg, source=ctx.providers.OpenStreetMap.Mapnik)
    ax.add_artist(ScaleBar(1, location="lower right"))

    north_arrow(ax=ax)

    return fig, ax


def north_arrow(ax: Axes, x: float = 0.95, y: float = 0.2, arrow_length: float = 0.1) -> Axes:
    """
    Adds a north arrow to a plot in axis ax

    Args:
        ax (Axes): axes to plot north arrow on
        x (float): relative location of arrow on x-axis (0-1).
        y (float): relative location of arrow on y-axis (0-1).
        arrow_length (float): size of arrow.

    Returns:
        ax (Axes): axes on which north arrow was plotted
    """
    ax.annotate(
        "N",
        xy=(x, y),
        xytext=(x, y - arrow_length),
        arrowprops=dict(facecolor="black", width=5, headwidth=15),
        ha="center",
        va="center",
        fontsize=20,
        xycoords=ax.transAxes,
    )
    return ax
