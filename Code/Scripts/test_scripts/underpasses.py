import sys
from copy import copy

import geopandas as gpd
import networkx as nx
import numpy as np
from shapely import affinity
from shapely.geometry import LineString, Point

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

        # out_gdf.loc[name, "h_date"] = -10

    return out_gdf


if __name__ == "__main__":
    ahn_path = r"D:\Work\Project\P1414\GIS\AHN\AHN4_WSS_filled.tif"
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
    wl_gdf = wl_gdf.loc[~wl_gdf["tunnel"], :]
    kering_gdf = gpd.read_file(r"D:\Work\Project\P1414\GIS\HYDAMO\Combined_keringen.gpkg")
    tunnel_gdf = gpd.read_file(r"D:\Work\Project\P1414\GIS\Keringen_met_hoogte\tunnel.shp")

    rd_gdf["fysiekvoor"] = rd_gdf["fysiekvoor"].fillna("")
    brug_bool = rd_gdf["fysiekvoor"].str.contains("brug")
    tunnel_bool = rd_gdf["fysiekvoor"].str.contains("tunnel")
    knoop_bool = rd_gdf["fysiekvoor"].str.contains("knooppunt")

    comb_bool = brug_bool | tunnel_bool  # | knoop_bool
    rd_gdf = rd_gdf.loc[~comb_bool, :]

    gdf["fysiekvoor"] = gdf["fysiekvoor"].fillna("")
    brug_bool = gdf["fysiekvoor"].str.contains("brug")
    tunnel_bool = gdf["fysiekvoor"].str.contains("tunnel")
    knoop_bool = gdf["fysiekvoor"].str.contains("knooppunt")

    comb_bool = brug_bool  # & ~knoop_bool

    in_gdf = gdf.loc[comb_bool, :]
    G = nx.Graph(gdf_to_nx(gdf_network=in_gdf))
    H = combine_straight_branches(G=G)
    out_branches_list = []
    for x, y, data in H.edges(data=True):
        out_branches_list.append(data)
    brug_gdf = gpd.GeoDataFrame(data=out_branches_list, geometry="geometry", crs=in_gdf.crs)

    brug_gdf = brug_gdf.loc[
        (brug_gdf["geometry"].length < 500) & (brug_gdf["geometry"].length > 10), :
    ]

    wl_gdf["geometry"] = wl_gdf["geometry"].buffer(25, cap_style=2)
    cols = brug_gdf.columns
    brug_gdf = brug_gdf.sjoin(wl_gdf, how="left")
    brug_gdf = brug_gdf.loc[brug_gdf["code"].isna(), cols]
    brug_gdf.to_file(
        r"D:\Work\Project\P1414\GIS\Wegen\Top10NLwegen_en_spoorwegen_clipped_merged_brug.shp"
    )

    rd_out_gdf = rd_gdf.sjoin(brug_gdf, how="left", predicate="crosses")
    rd_out_gdf = rd_out_gdf.loc[~rd_out_gdf["id_right"].isna(), :]

    rd_out_gdf["geometry"] = rd_out_gdf["geometry"].buffer(100)
    rd_out_gdf = rd_out_gdf.dissolve(by="hoogtenive_left").explode()

    rd_out_gdf.to_file(r"D:\Work\Project\P1414\GIS\Wegen\underpass_buffers.shp")

    blen = np.sqrt(2 * 500**2) + 1
    interp_range = 0.1
    in_gdf = copy(brug_gdf)
    u_list = []
    for name, bridge in in_gdf.iterrows():
        if name not in brug_gdf.index:
            print("skipped " + str(name))
            continue

        midpoint = bridge.geometry.interpolate(0.5, normalized=True)
        rotated_bridge = affinity.rotate(bridge.geometry, angle=90, origin=midpoint)

        llen = rotated_bridge.length
        rescale_factor = blen / llen

        p1 = rotated_bridge.interpolate(0.5 - interp_range / 2, normalized=True)
        p2 = rotated_bridge.interpolate(0.5 + interp_range / 2, normalized=True)

        if p2.x > p1.x:
            dx_o = p2.x - p1.x
            dy_o = p2.y - p1.y
        else:
            dx_o = p1.x - p2.x
            dy_o = p1.y - p2.y

        if dx_o != 0:
            s = dy_o / dx_o
            # given a line slope dy/dx (i.e. a/b), pythagoras a**2 + b**2 = c**2, and a desired length c
            # we can determine that b_n = c_n/sqrt(1+(a/b)**2) and a_n = b_n * (a/b)
            dx_1 = blen / np.sqrt(1 + s**2)
            dy_1 = dx_1 * s

        else:
            dx_1 = 0
            dy_1 = blen

        # centroid = rotated_line.centroid
        c_x, c_y = midpoint.x, midpoint.y

        offset_l = offset_r = 0.5

        list_of_points = [
            Point(c_x - offset_l * dx_1, c_y - offset_l * dy_1),
            Point(c_x + offset_r * dx_1, c_y + offset_r * dy_1),
        ]

        n_bridge = copy(bridge)
        # n_bridge["geometry"] = rotated_bridge
        n_bridge["geometry"] = LineString(list_of_points)
        bool_crosses = brug_gdf["geometry"].crosses(n_bridge["geometry"])
        brug_gdf = brug_gdf.loc[~bool_crosses, :]
        u_list.append(n_bridge.to_dict())

    u_gdf = gpd.GeoDataFrame(data=u_list, geometry="geometry", crs=in_gdf.crs)
    u_gdf.to_file(r"D:\Work\Project\P1414\GIS\Wegen\Underpasses.shp")

    kering_gdf["geometry"] = kering_gdf["geometry"].buffer(2.5)
    u_gdf = u_gdf.overlay(kering_gdf, how="difference").explode()

    tunnel_gdf["geometry"] = tunnel_gdf["geometry"].buffer(5)
    u_gdf = u_gdf.overlay(tunnel_gdf, how="difference").explode()

    u_gdf.to_file(r"D:\Work\Project\P1414\GIS\Wegen\Underpasses_overlayed.shp")
    out_u_gdf = copy(u_gdf)
    for name, upass in u_gdf.iterrows():
        bool_cross = in_gdf["geometry"].crosses(upass["geometry"])
        if np.sum(bool_cross.values[:]) == 0:
            out_u_gdf.drop(index=name, inplace=True)
            print("dropped " + str(name))

    out_u_gdf = add_height_to_linestrings(gdf=out_u_gdf, ahn_path=ahn_path, buffer=11)
    out_u_gdf = add_tunnel_dims(gdf=out_u_gdf)
    out_u_gdf = add_height_to_linestrings(gdf=out_u_gdf, ahn_path=ahn_path, buffer=10)
    out_u_gdf.to_file(r"D:\Work\Project\P1414\GIS\Wegen\Underpasses_filtered_500m.shp")
