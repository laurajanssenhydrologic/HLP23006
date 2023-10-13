import geopandas as gpd
import networkx as nx
import numpy as np

from add_ahn_height_to_fw import add_height_to_linestrings
from network_tools import combine_straight_branches, gdf_to_nx

if __name__ == "__main__":
    ahn_path = r"D:\Work\Project\P1414\GIS\AHN\AHN_merged.tif"
    input_path = r"D:\Work\Project\P1414\GIS\TOP10NLRelief\relief_clipped.shp"
    buffer_path = r"D:\Work\Project\P1414\GIS\HYDAMO\Combined_test_v7.3.gpkg"
    output_path = r"D:\Work\Project\P1414\GIS\Keringen_met_hoogte\relief_v2.shp"

    gdf = gpd.read_file(input_path).to_crs(crs="EPSG:28992").explode()
    b_gdf = gpd.read_file(buffer_path, layer="keringen").to_crs(crs="EPSG:28992")
    b_gdf["geometry"] = b_gdf["geometry"].buffer(15)

    in_gdf = gpd.overlay(gdf, b_gdf, how="difference").explode()

    G = nx.Graph(gdf_to_nx(gdf_network=in_gdf))
    H = combine_straight_branches(G=G)

    out_branches_list = []
    for x, y, data in H.edges(data=True):
        out_branches_list.append(data)

    _gdf = gpd.GeoDataFrame(data=out_branches_list, geometry="geometry", crs=in_gdf.crs)

    out_gdf = add_height_to_linestrings(gdf=_gdf, ahn_path=ahn_path, buffer=11)
    out_gdf.to_file(output_path)
