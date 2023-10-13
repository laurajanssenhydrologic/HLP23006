import sys
from copy import copy

import geopandas as gpd
import networkx as nx
import numpy as np

sys.path.append("D:\Work\git\GIS_tools\Code")
from geo_tools.add_ahn_height_to_fw import add_height_to_linestrings
from geo_tools.network_tools import combine_straight_branches
from geo_tools.networkx_tools import gdf_to_nx


def add_tunnel_dims(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    rail_width = [10, 15, 20]
    road_width = 3.5
    gdf["aantalspor"] = gdf["aantalspor"].str.strip()
    out_gdf = copy(gdf)
    for name, row in gdf.iterrows():
        if row["aantalspor"] == "enkel":
            out_gdf.loc[name, "width"] = rail_width[0]
            out_gdf.loc[name, "height"] = 4.7

        elif row["aantalspor"] == "dubbel":
            out_gdf.loc[name, "width"] = rail_width[1]
            out_gdf.loc[name, "height"] = 4.7

        elif row["aantalspor"] == "meervoudig":
            out_gdf.loc[name, "width"] = rail_width[2]
            out_gdf.loc[name, "height"] = 4.7

        else:
            if not np.isnan(row["aantalrijs"]):
                out_gdf.loc[name, "width"] = row["aantalrijs"] * road_width

            else:
                if row["verharding"] == "> 7 meter":
                    out_gdf.loc[name, "width"] = 7
                elif row["verharding"] == "4 - 7 meter":
                    out_gdf.loc[name, "width"] = (4 + 7) / 2

            out_gdf.loc[name, "height"] = 4.2

        out_gdf.loc[name, "h_date"] = -10

    return out_gdf


if __name__ == "__main__":
    ahn_path = r"D:\Work\Project\P1414\GIS\AHN\AHN_merged.tif"
    brug_output_path = r"D:\Work\Project\P1414\GIS\Keringen_met_hoogte\brug.shp"
    input_path = r"D:\Work\Project\P1414\GIS\Wegen\Top10NLwegen_en_spoorwegen_clipped.shp"
    output_path = (
        r"D:\Work\Project\P1414\GIS\Wegen\Top10NLwegen_en_spoorwegen_clipped_merged_brug.shp"
    )
    output_tunnel_path = r"D:\Work\Project\P1414\GIS\Keringen_met_hoogte\tunnel.shp"

    gdf = gpd.read_file(input_path).to_crs(crs="EPSG:28992").explode()
    rd_gdf = gpd.read_file("D:\Work\Project\P1414\GIS\Wegen\TOP10NLwegen_clipped.shp")
    wl_gdf = gpd.read_file(
        "D:\Work\Project\P1414\GIS\HYDAMO\Combined_test_v14_WBD.gpkg", layer="waterloop"
    ).to_crs(crs="EPSG:28992")

    rd_gdf["fysiekvoor"] = rd_gdf["fysiekvoor"].fillna("")
    brug_bool = rd_gdf["fysiekvoor"].str.contains("brug")
    tunnel_bool = rd_gdf["fysiekvoor"].str.contains("tunnel")
    knoop_bool = rd_gdf["fysiekvoor"].str.contains("knooppunt")

    comb_bool = brug_bool | tunnel_bool | knoop_bool
    rd_gdf = rd_gdf.loc[~comb_bool, :]

    gdf["fysiekvoor"] = gdf["fysiekvoor"].fillna("")
    brug_bool = gdf["fysiekvoor"].str.contains("brug")
    tunnel_bool = gdf["fysiekvoor"].str.contains("tunnel")
    knoop_bool = gdf["fysiekvoor"].str.contains("knooppunt")

    comb_bool = brug_bool | ~knoop_bool

    print(np.sum(comb_bool), np.sum(brug_bool), np.sum(tunnel_bool))

    in_gdf = gdf.loc[comb_bool, :]
    G = nx.Graph(gdf_to_nx(gdf_network=in_gdf))
    H = combine_straight_branches(G=G)
    out_branches_list = []
    for x, y, data in H.edges(data=True):
        out_branches_list.append(data)
    out_gdf = gpd.GeoDataFrame(data=out_branches_list, geometry="geometry", crs=in_gdf.crs)

    out_gdf = out_gdf.loc[
        (out_gdf["geometry"].length < 500) & (out_gdf["geometry"].length > 10), :
    ]

    # out_gdf["geometry"] = out_gdf["geometry"].buffer(25, cap_style=2)
    # out_gdf = out_gdf.dissolve(by=None).explode()

    wl_gdf["geometry"] = wl_gdf["geometry"].buffer(5, cap_style=2)
    cols = out_gdf.columns
    out_gdf = out_gdf.sjoin(wl_gdf, how="left")
    out_gdf = out_gdf.loc[out_gdf["code"].isna(), cols]
    out_gdf.to_file(
        r"D:\Work\Project\P1414\GIS\Wegen\Top10NLwegen_en_spoorwegen_clipped_merged_brug.shp"
    )
    rd_out_gdf = rd_gdf.sjoin(out_gdf, how="left", predicate="crosses")
    rd_out_gdf = rd_out_gdf.loc[~rd_out_gdf["id_right"].isna(), :]

    rd_out_gdf["geometry"] = rd_out_gdf["geometry"].buffer(25)
    rd_out_gdf = rd_out_gdf.dissolve(by="hoogtenive_left").explode()

    rd_out_gdf.to_file(
        r"D:\Work\Project\P1414\GIS\Wegen\Top10NLwegen_en_spoorwegen_clipped_merged_underpasses.shp"
    )


    
    # in_gdf = gpd.GeoDataFrame(data=out_branches_list, geometry="geometry", crs=in_gdf.crs)

    # out_gdf = add_height_to_linestrings(gdf=in_gdf, ahn_path=ahn_path, buffer=11)
    # out_gdf.to_file(output_path)

    # brug_gdf = gdf.loc[brug_bool, :]
    # out_brug_gdf = add_height_to_linestrings(gdf=brug_gdf, ahn_path=ahn_path)
    # out_brug_gdf = add_tunnel_dims(gdf=out_brug_gdf)
    # out_brug_gdf.to_file(brug_output_path)

    # tunnel_gdf = gdf.loc[tunnel_bool, :]
    # out_tunnel_gdf = add_height_to_linestrings(gdf=tunnel_gdf, ahn_path=ahn_path)
    # out_tunnel_gdf = add_tunnel_dims(gdf=out_tunnel_gdf)
    # out_tunnel_gdf.to_file(output_tunnel_path)
