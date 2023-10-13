import csv
from typing import List

import numpy as np
import pandas as pd
import rasterio
import xarray as xr
from rasterio.crs import CRS
from rasterio.errors import CRSError
from rasterio.transform import AffineTransformer, from_origin
from rasterio.windows import from_bounds
from tqdm import tqdm


def write_tiff(
    output_file_path: str, new_grid_data: np.ndarray, transform: AffineTransformer, epsg: int
) -> None:
    """
    Saves new_grid_data to a tiff file. new_grid_data should be a raster.
    Based on HydroLogic Inundation Toolbox, part of HYDROLIB.

    Args:
        output_file_path (str): location of the output file
        new_grid_data (np.ndarray): data at grid points
        transform (AffineTransformer): transform of the raster
        epsg (int): coordinate reference system (CPS) that is stored in the tiff-file

    Returns:
        None
    """

    # obtain target crs, if not given use EPSG:28992 (Dutch RDS)
    try:
        raster_crs = CRS.from_epsg(epsg)
    except CRSError:
        raster_crs = CRS.from_epsg(28992)

    # Create tif file with relevant options
    t_file = rasterio.open(
        output_file_path,
        "w",
        driver="GTiff",
        height=new_grid_data.shape[0],
        width=new_grid_data.shape[1],
        count=1,
        dtype=new_grid_data.dtype,
        crs=raster_crs,
        transform=transform,
    )

    # Write and close file
    t_file.write(new_grid_data, 1)
    t_file.close()


def raster_transform(
    source: rasterio.io.DatasetReader,
    destination: rasterio.io.DatasetReader,
    output_path: str = None,
    resampling: int = 5,
    epsg: int = 28992,
    nik_to_man: bool = False,
) -> None:
    """
    Transforms source raster shape to destination raster shape.
    Output is saved at output_path. Resampling options from rasterio are available.

    Args:
        source (rasterio.io.DatasetReader): source raster from which data is extracted.
        destination (asterio.io.DatasetReader): destination raster to which's shape the new raster will conform.
        output_path (str): location to save the output raster. Default None will result in no files being saved to disk
        resampling (int): rasterio resampling options (default: 5 = average). Options are:
                nearest = 0,
                bilinear = 1,
                cubic = 2,
                cubic_spline = 3,
                lanczos = 4,
                average = 5,
                mode = 6,
                gauss = 7,
                max = 8,
                min = 9,
                med = 10,
                q1 = 11,
                q3 = 12,
                sum = 13,
                rms = 14.
        epsg (int): coordinate reference system. Default is Dutch RDS.
        nik_to_man (bool): Whether to convert roughness from Nikuradse to Manning. Defaults to False


    Returns:
        s_data_np (np.ndarray): transformed raster in numpy
        s_data_da (xr.DataArray): transformed raster in Xarray DataArray format
    """
    # Read Rasterio transforms from both source and destination rasters
    s_transform = source.transform
    d_transform = destination.transform

    # Read destination data for its shape, and create a new source data array with the same shape
    d_data = destination.read(1)
    s_data = np.zeros(d_data.shape)

    # Obtain coordinate boundary of the pixel (upper left and lower right bounds)
    ul_x, ul_y = rasterio.transform.xy(transform=d_transform, rows=0, cols=0, offset="ul")
    lr_x, lr_y = rasterio.transform.xy(
        transform=d_transform, rows=d_data.shape[0], cols=d_data.shape[1], offset="lr"
    )

    # Create a Rasterio window from pixel bounds and source_transform
    window = from_bounds(left=ul_x, top=ul_y, right=lr_x, bottom=lr_y, transform=s_transform)

    # Read source data within window and resample
    s_data_np = source.read(
        out_shape=(1, d_data.shape[0], d_data.shape[1]), window=window, resampling=resampling
    )[0, :, :]

    if nik_to_man:
        s_data_np = roughness_nikuradse_to_manning(s_data_np)

    xs, _ = rasterio.transform.xy(
        transform=d_transform, rows=0, cols=np.arange(d_data.shape[1]), offset="center"
    )
    _, ys = rasterio.transform.xy(
        transform=d_transform, rows=np.arange(d_data.shape[0]), cols=0, offset="center"
    )

    # store data in DataArray
    s_data_da = xr.DataArray(
        data=s_data_np,
        dims=["y", "x"],
        coords={"x": xs, "y": ys},
    )

    # If output path is given, save data to disk
    if output_path is not None:
        if (".tif" in output_path.lower()) or (".tiff" in output_path.lower()):
            # write tiff
            write_tiff(
                output_file_path=output_path,
                new_grid_data=s_data_np,
                transform=d_transform,
                epsg=epsg,
            )
        elif ".nc" in output_path.lower():
            s_data_da.to_netcdf(path=output_path, mode="w", format="NETCDF4")
        else:
            raise NameError("File not supported")

    return s_data_np, s_data_da


def raster_to_xyz(raster_path, xyz_output_path):
    with rasterio.open(raster_path) as src:
        # read image
        image = src.read()
        # transform image
        bands, rows, cols = np.shape(image)
        image1 = image.reshape(rows * cols, bands)

        # bounding box of image
        l, b, r, t = src.bounds
        # resolution of image
        res = src.res

        # meshgrid of X and Y
        x = np.arange(l, r, res[0])
        y = np.arange(t, b, -res[0])
        X, Y = np.meshgrid(x, y)

        # flatten X and Y
        newX = X.flatten()
        newY = Y.flatten()

        # join XY and Z information
        export = np.column_stack((newX, newY, image1))

        with open(xyz_output_path, "w", newline="") as fp:
            a = csv.writer(fp, delimiter=" ")
            a.writerows(export)
            fp.close()  # close file


def roughness_nikuradse_to_manning(raster):
    manning_raster = (raster ** (1 / 6)) * (1 / 25)

    manning_raster[raster == 0] = 0.01  # glass
    manning_raster = np.nan_to_num(manning_raster, nan=0.023)  # default value
    return manning_raster


def create_roughness_xyz(
    xnodes,
    ynodes,
    dx,
    dy,
    roughness_tif_path,
    roughness_mesh_name,
    roughness_xyz_name,
    convert_to_manning=False,
):

    gridX, _ = np.meshgrid(np.unique(xnodes), np.unique(ynodes))

    grid_data = np.zeros_like(gridX)

    affine_transform = from_origin(np.min(xnodes), np.max(ynodes), dx, dy)

    # Write a tif with the 2D meshgrid
    write_tiff(r"Testfile_str.tif", grid_data, affine_transform, 28992)

    # Open the newly created mesh tif and the roughness tif and resample
    grid_tif = rasterio.open(r"Testfile_str.tif")
    roughness = rasterio.open(roughness_tif_path)

    # _, _ = raster_transform(roughness, grid_tif, roughness_mesh_name)
    _, _ = raster_transform(
        roughness, grid_tif, roughness_mesh_name, nik_to_man=convert_to_manning
    )

    raster_to_xyz(roughness_mesh_name, roughness_xyz_name)


# raster_to_xyz(r"D:\work\P1414_ROI\GIS\roughness_meshgrid.tif", r"D:\work\P1414_ROI\GIS\roughness_list.xyz")
