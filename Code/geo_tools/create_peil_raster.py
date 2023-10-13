from copy import copy

import geopandas as gpd
import numpy as np
import rasterio
from rasterio import features

boezem_gdf = gpd.read_file(r"D:\Work\Project\P1414\GIS\peilen\boezem_v2.shp")
orig_raster = r"D:\Work\Project\P1414\GIS\peilen\wd_0_v3.tif"
outp_raster = r"D:\Work\Project\P1414\GIS\peilen\wd_0_v4.tif"


src = rasterio.open(orig_raster)
meta = src.meta.copy()
meta.update(compress="lzw")

with rasterio.open(outp_raster, "w+", **meta) as out:
    out_arr = out.read(1)

    # this is where we create a generator of geom, value pairs to use in rasterizing
    shapes = ((geom, value) for geom, value in zip(boezem_gdf.geometry, boezem_gdf.peil))

    burned = features.rasterize(shapes=shapes, fill=np.nan, out=out_arr, transform=out.transform)
    burned[burned < -2] = -10
    out.write_band(1, burned)
