import os
import sys

import matplotlib.pyplot as plt
import numpy as np

currentdir = os.path.dirname(os.getcwd())
sys.path.append(r"D:\Work\git\GIS_tools\HydroLogic_Inundation_Toolbox\Readers")

from flowmeshreader import load_fou_data, mesh_to_tiff
from plotting import raster_plot_with_context

# set paths
model_path = r"D:\Work\Project\P1414\Models\Combined\V23\V23_WBD_500_exported\dflowfm\output"
input_file_path = model_path + r"\\DFM_fou.nc"
output_file_path = model_path + r"\\max_depth.tiff"
output_png_file_path = model_path + r"\\max_depth.png"

# raster options
resolution = 50  # m
dhydro_resolution = 500  # m
distance_tol = np.ceil(np.sqrt(2 * dhydro_resolution**2))  # m
interpolation = r"nearest"

variable = r"Mesh2d_fourier001_max_depth"

# load mesh coordinates and data from netCDF
node_data = load_fou_data(input_file_path, variable)
node_data[node_data < 0.02] = np.nan

# convert to raster and save as tiff
_, _, grid_data = mesh_to_tiff(
    node_data,
    input_file_path,
    output_file_path,
    resolution,
    distance_tol,
    interpolation=interpolation,
)

fig, ax = raster_plot_with_context(
    raster_path=output_file_path,
    epsg=28992,
    clabel="water depth (m)",
    cmap="Reds",
    title="Maximum water depth",
)

fig.savefig(output_png_file_path, dpi=300)
