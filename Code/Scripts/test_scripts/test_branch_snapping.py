from copy import copy

import geopandas as gpd
import matplotlib.pyplot as plt

folder = r"D:\Work\Project\P1414"
gpkg_path = folder + r"\GIS\HYDAMO\combined_test.gpkg"

branches_gdf = gpd.read_file(gpkg_path, layer="waterloop")
selection_list = [
    "56c8f45f-a109-4606-b61d-bb164412b914",
    "3a1eb578-c8f5-4c78-949e-ac5611fcc9c8",
    "466913ac-8557-4c03-9c1d-13c8b3a81d8b",
]
selected_branches = branches_gdf[branches_gdf["globalid"].isin(selection_list)]
print(selected_branches)


branch1 = branches_gdf[branches_gdf["globalid"] == "56c8f45f-a109-4606-b61d-bb164412b914"]
branch2 = branches_gdf[branches_gdf["globalid"] == "3a1eb578-c8f5-4c78-949e-ac5611fcc9c8"]
branch3 = branches_gdf[branches_gdf["globalid"] == "466913ac-8557-4c03-9c1d-13c8b3a81d8b"]
union_result = branch1.geometry.iloc[0].union(branch2.geometry.iloc[0])
union_result = selected_branches.geometry.unary_union
union_result = branches_gdf.geometry.unary_union

# print(union_result)
buffered_branches = copy(branches_gdf)
buffered_branches.geometry = buffered_branches.geometry.buffer(1)
gdf = gpd.GeoDataFrame(union_result, columns=["geometry"], geometry="geometry", crs=28992)
intersected = gdf.sjoin(buffered_branches, how="left", predicate="within")
# uniques = gpd.GeoDataFrame(
#     intersected.drop_duplicates(subset="geometry"), geometry="geometry", crs=28992
# )
intersected.plot(column="globalid", cmap="Set1")
intersected.to_file(r"D:\Work\Project\P1414\GIS\test.shp")
plt.show()
print(branches_gdf.shape)
print(intersected.shape)

print(branch1.touches(branch2, align=False))
print(branch3.touches(branch2, align=False))


branch1_points = list(branch1.geometry.iloc[0].coords)
branch2_points = list(branch2.geometry.iloc[0].coords)
branch3_points = list(branch3.geometry.iloc[0].coords)

branch1_points = [branch1_points[0], branch1_points[-1]]
branch2_points = [branch2_points[0], branch2_points[-1]]
branch3_points = [branch3_points[0], branch3_points[-1]]

if (branch1_points[0] in branch2_points) or (branch1_points[-1] in branch2_points):
    print("branch 1 en 2 zitten aan elkaar")
if (branch3_points[0] in branch2_points) or (branch3_points[-1] in branch2_points):
    print("branch 2 en 3 zitten aan elkaar")

ip1 = branch1.geometry.iloc[0].intersection(branch2.geometry.iloc[0])
ip2 = branch3.geometry.iloc[0].intersection(branch2.geometry.iloc[0])
ip3 = branch3.geometry.iloc[0].intersection(branch1.geometry.iloc[0])
print(ip1)
print(ip2)
print(ip3)
print(ip3.is_empty)
