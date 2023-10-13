import numpy as np
import rasterio
from rasterio.crs import CRS

input_path = r"D:\Work\Project\P1414\GIS\AHN\AHN4_WSS.tif"
output_path = r"D:\Work\Project\P1414\GIS\AHN\AHN4_WSS_filled.tif"

raster = rasterio.open(input_path)
r_array = raster.read(1)
print(np.sum(np.isnan(r_array)))
r_array[r_array > 30] = -10

r_file = rasterio.open(
    output_path,
    "w",
    driver="GTiff",
    height=r_array.shape[0],
    width=r_array.shape[1],
    count=1,
    dtype=r_array.dtype,
    crs=CRS.from_epsg(28992),
    transform=raster.transform,
)
# Write and close file
r_file.write(r_array, 1)
r_file.close()
