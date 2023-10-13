import geopandas as gpd
import matplotlib.pyplot as plt

# gpkg_path = "D:\Work\Project\P1414\GIS\HYDAMO\combined_test.gpkg"
folder = r"D:\Work\Project\P1414"
gpkg_path = folder + r"\GIS\HYDAMO\HDSR_selection.gpkg"

branches_gdf = gpd.read_file(gpkg_path, layer="waterloop")
culvert_gdf = gpd.read_file(gpkg_path, layer="duiker")


def checks(branches_gdf, culvert):
    _branch = branches_gdf.clip(
        gpd.GeoSeries(culvert.buffer(0.1, cap_style=2), crs=branches_gdf.crs),
        keep_geom_type=True,
    )
    if _branch.shape[0] > 0:
        branch = _branch.geometry.values[0].buffer(0.1)

        contains = culvert.covered_by(branch)
        crosses = branch.crosses(culvert)
        overlaps = branch.overlaps(culvert)
        touches = branch.touches(culvert)
        within = culvert.within(branch)
        print(contains, crosses, overlaps, touches, within)


# culvert = culvert_gdf[culvert_gdf["code"] == "SY6012"].geometry.values[0]

# checks(branches_gdf, culvert)


# culvert = culvert_gdf[culvert_gdf["code"] == "D3428"].geometry.values[0]
# checks(branches_gdf, culvert)

# for ix, culvert in tqdm(culvert_gdf.iterrows(), total=culvert_gdf.shape[0]):
#     checks(branches_gdf, culvert.geometry)

correct_culverts = culvert_gdf.overlay(
    gpd.GeoDataFrame(
        branches_gdf.dissolve(by=None).buffer(0.1),
        columns=["geometry"],
        crs=branches_gdf.crs,
        geometry="geometry",
    ),
    how="intersection",
)
print(correct_culverts)

base = branches_gdf.plot()
correct_culverts.plot(ax=base, color="red")
plt.show()
