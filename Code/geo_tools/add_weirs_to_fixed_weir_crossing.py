import geopandas as gpd
import numpy as np
from scipy.spatial import KDTree
from shapely.geometry import Point

gpkg_path = r"D:\Work\Project\P1414\GIS\HYDAMO\Combined_test_v15_Markermeer.gpkg"
max_distance = 10

branch_gdf = gpd.read_file(gpkg_path, layer="waterloop")
fw_gdf = gpd.read_file(gpkg_path, layer="keringen")
weirs_gdf = gpd.read_file(gpkg_path, layer="stuw")

bool_array = branch_gdf["code"].str.contains("tunn") | branch_gdf["code"].str.contains("unde")
branch_gdf = branch_gdf.loc[~bool_array, :]
branch_gdf["geometry"] = branch_gdf["geometry"].buffer(0.1)
fw_gdf["geometry"] = fw_gdf["geometry"].buffer(0.1)

new_weirs_gdf = gpd.overlay(branch_gdf, fw_gdf, how="intersection").explode()

print(new_weirs_gdf.shape)
print(new_weirs_gdf.head())

point_list_base = []
for ix, weir in weirs_gdf.iterrows():
    if isinstance(weir.geometry, Point):
        coords = weir.geometry.coords[0]
    else:
        coords = weir.geometry.centroid.coords[0]

    if len(coords) > 2:
        coords = coords[:2]

    point_list_base.append(coords)


point_list_match = []
for ix, weir in new_weirs_gdf.iterrows():
    if isinstance(weir.geometry, Point):
        coords = weir.geometry.coords[0]
    else:
        coords = weir.geometry.centroid.coords[0]

    if len(coords) > 2:
        coords = coords[:2]

    point_list_match.append(coords)

# build spatial KDTree from start and end points
kdtree_base = KDTree(point_list_base)
kdtree_match = KDTree(point_list_match)

# Set up a sparse distance matrix
sdm = kdtree_base.sparse_distance_matrix(
    kdtree_match, max_distance=max_distance, output_type="coo_matrix"
)

# Check for non-zero (explicit zeroes allowed) in base dataset
nnz_base = sdm.getnnz(axis=0)
bool_array = nnz_base > 0

print(np.sum(bool_array))
print(bool_array.shape)
print(new_weirs_gdf.shape)

new_weirs_gdf = new_weirs_gdf.loc[~bool_array, :]
new_weirs_gdf["geometry"] = new_weirs_gdf["geometry"].centroid
new_weirs_gdf.to_file(r"D:\Work\Project\P1414\GIS\Keringen\nieuwe_stuwen.shp")
