import geopandas as gpd
import numpy as np
import rasterio
from rasterio.features import shapes
from tqdm import tqdm

raster_path = r"D:\Work\Project\P1414\Models\Combined\V21_WBD_HD\run_model_changed_coords\dflowfm\output\fig_wl\480.tiff"
shape_path = r"D:\Work\Project\P1414\Models\Combined\V21_WBD_HD\run_model_changed_coords\dflowfm\output\fig_wl"


with rasterio.open(raster_path) as src:
    data = src.read(1).astype("float32")

    shp_iter = shapes(data, mask=None, transform=src.transform)
    pol_list = []
    for ix, (pol, val) in tqdm(enumerate(shp_iter)):
        if np.isnan(val):
            continue
        vals = dict([("maxlevel", val), ("minlevel", val), ("name", ix), ("step", 1)])
        pol = dict([("geometry", pol), ("properties", vals)])
        pol_list.append(pol)

gdf = gpd.GeoDataFrame.from_features(pol_list)
gdf["areas"] = gdf.geometry.area / 1e6
areas = gdf["areas"].to_numpy()
area = gdf.geometry.area.sort_values().sum() / 1e6
print(area)

max_area = 500
parts = np.ceil(area / max_area).astype(int)

start_ix = 0
for part in range(parts):
    gdf_path = shape_path + r"\shape_{}".format(part)
    _areas = areas[start_ix:]
    cum_areas = np.cumsum(_areas)
    bool_array = cum_areas < max_area
    stop_ix = np.sum(bool_array) + start_ix

    _gdf = gdf.iloc[start_ix:stop_ix, :]
    print(_gdf.areas.sum())

    _gdf = _gdf.drop(columns="areas")
    _gdf.to_file(gdf_path)
    start_ix = stop_ix


# gdf.to_file(
#     r"D:\Work\Project\P1414\Models\Combined\V21_WBD_HD\run_model_changed_coords\dflowfm\output\fig_wl\480.shp"
# )
