from copy import copy

import geopandas as gpd
import numpy as np
import rasterio
from shapely.geometry import LineString

input_path = r"D:\Work\Project\P1414\GIS\Keringen_met_hoogte\hdsr.shp"


gdf = gpd.read_file(input_path).to_crs(crs="EPSG:28992").explode()
coord_list = gdf.iloc[0, :].geometry.coords[:]
print(coord_list[0])
