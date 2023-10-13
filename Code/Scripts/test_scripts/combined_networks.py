from itertools import compress

import geopandas as gpd
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from scipy.spatial import KDTree
from shapely.geometry import LineString, Point

# hdsr_gdf = gpd.read_file(
#     # r"D:\Work\Project\P1414\GIS\HDSR\Legger\Hydro_Objecten(2)\HydroObject.shp"
#     r"D:\Work\Project\P1414\GIS\Uitgesneden watergangen\HDSR_v00_test.shp"
# ).explode(ignore_index=True)
# hdsr_gdf = hdsr_gdf.loc[hdsr_gdf["IWS_W_WATB"] > 5, :]
# # hdsr_gdf = hdsr_gdf.loc[hdsr_gdf["NAAM"] != "Vecht", :]
# hdsr_gdf = hdsr_gdf.loc[hdsr_gdf["NAAM"] != "Nederrijn", :]
# hdsr_gdf = hdsr_gdf.loc[hdsr_gdf["NAAM"] != "Lek", :]
# hdsr_gdf = hdsr_gdf.loc[hdsr_gdf["NAAM"] != "Amsterdam - Rijn Kanaal", :]
# print("HDSR done")
hhd_gdf = gpd.read_file(
    r"D:\Work\Project\P1414\GIS\Uitgesneden watergangen\HHD_v00_test.shp"
).explode(ignore_index=True)
hhr_gdf = gpd.read_file(
    r"D:\Work\Project\P1414\GIS\Uitgesneden watergangen\HHR_v00_test.shp"
).explode(ignore_index=True)
hhr_gdf = hhr_gdf.loc[hhr_gdf["BODEMBREED"] > 5, :]
# hhsk_gdf = gpd.read_file(
#     r"D:\Work\Project\P1414\GIS\Uitgesneden watergangen\HHSK_v00_test.shp"
# ).explode(ignore_index=True)
# hhsk_gdf = hhsk_gdf.loc[hhsk_gdf["WATERBREED"] > 5, :]

# rt_gdf = gpd.read_file(
#     r"D:\Work\Project\P1414\GIS\Rijntakken\Region_Network_Branches.shp"
# ).explode(ignore_index=True)

for ix, row in hhr_gdf.iterrows():
    if len(row.geometry.coords[0]) > 2:
        geometry = [Point(x, y) for x, y, _ in row.geometry.coords]
        hhr_gdf.loc[ix, "geometry"] = LineString(geometry)


# wagv_gdf = gpd.read_file(
#     r"D:\Work\Project\P1414\GIS\Uitgesneden watergangen\AGV_v00_test.shp"
# ).explode(ignore_index=True)
# wagv_gdf = wagv_gdf.loc[wagv_gdf["IWS_W_WATB"] > 5, :]
# print("WAGV done")


# buffered_hdsr = gpd.GeoDataFrame(
#     geometry=hdsr_gdf.dissolve(by=None).buffer(distance=25, cap_style=1)
# )

# # hdsr_gdf = hdsr_gdf.overlay(buffered_wagv, how="difference").explode()
# wagv_gdf = wagv_gdf.overlay(buffered_hdsr, how="difference").explode(keep_geom_type=True)
# print("HDSR fixed")

gdf1 = hhr_gdf
gdf2 = hhd_gdf
# skip branches further than max_dist away from another branch
gdf_list = [gdf1, gdf2]
kdtree_list = []
for gdf in gdf_list:
    point_list = []
    for ix, branch in gdf.iterrows():
        coords = branch.geometry.coords
        point_list.append(coords[0])
        point_list.append(coords[-1])

    # build spatial KDTree from start and end points
    _kdtree = KDTree(point_list)
    kdtree_list.append(_kdtree)

sdm = kdtree_list[0].sparse_distance_matrix(
    kdtree_list[1], max_distance=100, output_type="coo_matrix"
)

# matrix is symetrical, so picking one of the two axis to count number of non empty entries (explicit zeros allowed)
nnz = sdm.getnnz(axis=0)
print(*nnz[nnz != 0])
print(list(compress(point_list, nnz != 0)))

fig = plt.figure()
ax = fig.gca()
gdf1.plot(ax=ax, color="r")
gdf2.plot(ax=ax, color="g")

if np.sum(nnz != 0) > 0:
    x, y = zip(*list(compress(point_list, nnz != 0)))
    plt.scatter(x=x, y=y, color="b", s=10)
plt.show()
