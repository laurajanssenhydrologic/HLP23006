import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd

gpkg_path = r"D:\Work\Project\P1414\GIS\HYDAMO\Combined_test_v16_Markermeer.gpkg"
peil_list = [-0.43, -0.625, -0.4, -0.3, -0.45, 0.58, -1.05]
top10nl_water = r"D:\Work\GIS\TOP10NL\top10nl_Waterdeel.gpkg"

branch_gdf = gpd.read_file(gpkg_path, layer="waterloop")
fw_gdf = gpd.read_file(gpkg_path, layer="keringen")
peil_gdf = gpd.read_file(gpkg_path, layer="peilgebieden")
top10nlwater = gpd.read_file(top10nl_water, layer="top10nl_waterdeel_vlak")


branch_gdf = branch_gdf.loc[branch_gdf["peil"].isin(peil_list), :]
# branch_gdf["geometry"] = branch_gdf["geometry"].buffer(50)
branch_gdf = gpd.sjoin(branch_gdf, top10nlwater, how="right", predicate="intersects")
branch_gdf = branch_gdf.loc[~branch_gdf["peil"].isna(), :]

# branch_gdf_old = gpd.sjoin(branch_gdf, peil_gdf, how="left", predicate="intersects")
# branch_gdf = branch_gdf_old.loc[branch_gdf_old["peil_left"] == branch_gdf_old["peil_right"], :]
# branch_gdf = pd.concat(
#     [
#         branch_gdf,
#         branch_gdf_old.loc[
#             branch_gdf_old["code"].str.contains("wagv")
#             | branch_gdf_old["code"].str.contains("ark"),
#             :,
#         ],
#     ],
# )
# branch_gdf = branch_gdf.dissolve(by="peil").explode()


branch_gdf.to_file(r"D:\Work\Project\P1414\GIS\peilen\boezem_v2.shp")
branch_gdf.plot()


plt.show()
