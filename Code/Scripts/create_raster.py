#%%
import rasterio
from rasterio.transform import from_origin
import numpy as np

# Inladen van het raster
input_raster_path = 'D:\work\P23006\GIS\peilen\peilen_jp_25m_full.tif'
output_raster_path = 'D:\work\P23006\GIS\peilen\inipeilen_raster2d.tif'
#%%
with rasterio.open(input_raster_path) as src:
    # Initialiseer het raster met -10
    data = np.full((src.height, src.width), -10, dtype=src.dtypes[0]) 

    # Lees het raster en overschrijf de waarden die geen nodata zijn
    src_data = src.read(1)
    valid_data_mask = (src_data != src.nodata)
    #data[valid_data_mask] = src_data[valid_data_mask]

    # Verkrijg de transformatiematrix en pas deze aan als dat nodig is
    transform = src.transform
    # Als je de transformatiematrix wilt aanpassen, kun je dit hier doen

    # Oplslaan van het raster met de aangepaste waarden
    with rasterio.open(output_raster_path, 'w', driver='GTiff', height=src.height, width=src.width, count=1, dtype=str(data.dtype), crs=src.crs, transform=transform) as dst:
        dst.write(data, 1)
# %%
from shapely.geometry import Point
import rasterio
from rasterio.features import geometry_mask
import geopandas as gpd

name="129529.000000_462192.000000"
coordinates = name.split('_')
bound_coord = Point(float(coordinates[0]),float(coordinates[1]))

raster_path= r'D:\work\P23006\GIS\peilen\wd_0_v4_test_4_5m.tif'
with rasterio.open(raster_path) as src:
    # Create a mask for the point
    point_mask = geometry_mask([bound_coord], out_shape=src.shape, transform=src.transform, invert=True)

    # Read the raster values only where the point is located
    values_at_point = src.read(1, masked=True)[point_mask]
    waterlevel=values_at_point.data[0]


# %%
