import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from tqdm import tqdm

sys.path.append(r"D:\Work\git\GIS_tools\HydroLogic_Inundation_Toolbox\Readers")

from flowmeshreader import load_map_data, mesh_to_tiff
from plotting import raster_plot_with_context

# set paths
# output_folder = r"D:\Work\Project\P1414\Models\HHR\HM\run_fw_model\dflowfm\output"
output_folder = (
    r"D:\Work\Project\P1414\Models\Combined\V21_WBD_HD\run_model_changed_coords\dflowfm\output"
)
input_file_path = output_folder + r"\DFM_map.nc"


# raster options
resolution = 33  # m
dhydro_resolution = 100  # m
distance_tol = np.ceil(np.sqrt(2 * dhydro_resolution**2))  # m
interpolation = r"nearest"

for n in range(2):
    if n == 0:
        continue
        # # load mesh coordinates and data from netCDF
        variable = r"Mesh2d_waterdepth"
        var_str = "Water depth"
        cmap = "Reds"
        vmin = 0
        vmax = 10
        output_gif_path = output_folder + r"\wd.webp"
        output_fig_path = output_folder + r"\fig_wd"
        Path(output_fig_path).mkdir(exist_ok=True)
        map_data = load_map_data(input_file_path, variable)
        map_data[map_data < 0.01] = np.nan
    if n == 1:
        variable = r"Mesh2d_s1"
        var_str = "Water level"
        cmap = "RdBu_r"
        vmin = -10
        vmax = 10
        output_gif_path = output_folder + r"\wl.webp"
        output_fig_path = output_folder + r"\fig_wl"
        Path(output_fig_path).mkdir(exist_ok=True)
        map_data = load_map_data(input_file_path, variable)

    # Loop over frames and check if png exist.
    # If not, check if tiff exists to make png from.
    # If both don't exist, create both

    frames = []
    for ix in tqdm(range(480, 505, 25)):
        # for ix in tqdm(range(0, 72, 3)):
        # for ix in tqdm(range(1, 13, 1)):
        output_png_file_path = output_fig_path + r"\{}.png".format(ix)
        output_tiff_file_path = output_fig_path + r"\{}.tiff".format(ix)

        if Path(output_png_file_path).is_file():
            tqdm.write("{}.png found".format(ix))

        elif Path(output_tiff_file_path).is_file():
            tqdm.write("{}.tiff found".format(ix))
            fig, ax = raster_plot_with_context(
                raster_path=output_tiff_file_path,
                epsg=28992,
                clabel=var_str,
                cmap=cmap,
                title=var_str,
                vmin=vmin,
                vmax=vmax,
            )
            fig.savefig(output_png_file_path, dpi=300)
            plt.close(fig)

        else:
            tqdm.write("{}.tiff not found".format(ix))
            # convert to raster and save as tiff
            _, _, _ = mesh_to_tiff(
                map_data[ix, :],
                input_file_path,
                output_tiff_file_path,
                resolution,
                distance_tol,
                interpolation=interpolation,
            )

            fig, ax = raster_plot_with_context(
                raster_path=output_tiff_file_path,
                epsg=28992,
                clabel=var_str,
                cmap=cmap,
                title=var_str,
                vmin=vmin,
                vmax=vmax,
            )
            fig.savefig(output_png_file_path, dpi=300)
            plt.close(fig)

        new_frame = Image.open(output_png_file_path)
        frames.append(new_frame)

    frames[0].save(
        output_gif_path,
        format="webp",
        append_images=frames[1:],
        disposal=2,
        duration=500,
        loop=0,
        optimize=True,
        save_all=True,
        transparency=1,
    )
