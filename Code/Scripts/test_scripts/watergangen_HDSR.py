import os
import sys
import uuid
from copy import copy

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from shapely.geometry import LineString

sys.path.append(r"D:\Work\git\hl23006\Code\geo_tools")
from merge_networks import merge_networks
from network_tools import combine_straight_branches
from networkx_tools import gdf_to_nx


def combine_branches(in_gdf: gpd.GeoDataFrame, geometry_accuracy: float = 0) -> gpd.GeoDataFrame:
    in_gdf = in_gdf.explode(ignore_index=True)
    G = gdf_to_nx(in_gdf, decimals=geometry_accuracy)
    H = combine_straight_branches(G=G)
    out_branches_list = []
    for x, y, data in H.edges(data=True):
        out_branches_list.append(data)

    out_branches = gpd.GeoDataFrame(data=out_branches_list, geometry="geometry", crs=in_gdf.crs)
    return out_branches


def snap_endpoints(in_branches: gpd.GeoDataFrame, geometry_accuracy: float):
    in_branches["t_id"] = [str(uuid.uuid4()) for _ in range(in_branches.shape[0])]
    in_branches.set_index("t_id", inplace=True)
    out_branches = copy(in_branches)
    for ix, branch in in_branches.iterrows():
        point_list = []
        n_points = len(branch.geometry.coords[:])
        for jx, (x, y, *_) in enumerate(branch.geometry.coords[:]):
            if (jx == 0) or (jx == (n_points - 1)):
                point_list.append(
                    (np.around(x, geometry_accuracy), np.around(y, geometry_accuracy))
                )
            else:
                point_list.append((x, y))

        geometry = LineString(point_list)
        if geometry.length > 0:
            out_branches.loc[ix, "geometry"] = LineString(point_list)
        else:
            out_branches = out_branches.drop(index=ix)

    return out_branches


def snap_nodes(in_branches: gpd.GeoDataFrame, geometry_accuracy: float):
    """
    Snaps nodes to a defined accuracy

    Args:
        in_branches (gpd.GeoDataFrame): GeoDataFrame containing the branches
        geometry_accuracy (float): Number of decimals to round. This is done using numpy round.
        Negative values means rounding on the left side of the decimal marker.
        e.g. -1 would round to 10's, 0 to 1's, and 1 to 0.1's

    """
    in_branches["t_id"] = [str(uuid.uuid4()) for _ in range(in_branches.shape[0])]
    in_branches.set_index("t_id", inplace=True)
    out_branches = copy(in_branches)
    for ix, branch in in_branches.iterrows():

        point_list = []
        for x, y, *_ in branch.geometry.coords[:]:
            point_list.append((np.around(x, geometry_accuracy), np.around(y, geometry_accuracy)))

        geometry = LineString(point_list)
        if geometry.length > 0:
            out_branches.loc[ix, "geometry"] = LineString(point_list)
        else:
            out_branches = out_branches.drop(index=ix)

    return out_branches


def validate_network_topology(
    branches_gdf: gpd.GeoDataFrame, geometry_accuracy: float
) -> gpd.GeoDataFrame:
    """
    Function to validate the toplogy of the network.
    Specifically, this function ensures all connections between branches are valid
    and that each branch connects to another at the end of a branch.

    TODO: change to networkx

    Args:
        branches_gdf (gpd.GeoDataFrame): input geodataframe with branches (geometry should be linestrings)

    Returns:
        branches_gdf (gpd.GeoDataFrame): output geodataframe with branches that properly connect to one another
    """

    MLS_branches_bool = branches_gdf.geometry.type == "MultiLineString"
    if np.sum(MLS_branches_bool) > 0:
        print("warning: found multilinestrings")
        branches_gdf = branches_gdf.explode(ignore_index=True, index_parts=False)

    # First ensure that line-segments dont extend past connection points
    # this is done using shapely function unary_union
    union_result = branches_gdf.geometry.unary_union

    # Secondly, add the correct data from the origingal gdf to the new branches
    # this is done by adding data from old buffered branches to all new branches that fall within each polygon
    buffered_branches = copy(branches_gdf)
    buffered_branches.geometry = buffered_branches.geometry.buffer(0.01, cap_style=2)
    gdf = gpd.GeoDataFrame(union_result, columns=["geometry"], geometry="geometry", crs=28992)

    # add data from buffered branches if lines in gdf fall within. But keep geometry of branches in gdf
    intersected_gdf = gdf.sjoin(buffered_branches, how="left", predicate="within")

    # intersected_gdf = snap_nodes(in_branches=intersected_gdf, geometry_accuracy=geometry_accuracy)
    # if snap_distance is not None:
    #     if snap_distance > 0:
    #         intersected_gdf = delete_branches(gdf=intersected_gdf, snap_distance=snap_distance)
    #     else:
    #         raise ValueError("snap_distance should be > 0")

    return intersected_gdf


duikers_path = r"D:\Work\Project\HL-23006\GIS\Legger\Kokers_Lijnen.shp"
branches_path = (
    r"D:\Work\Project\HL-23006\GIS\Selectie_watergangen_21062023\Selectie_watergangen_21062023.shp"
)
np_path = r"D:\Work\Project\HL-23006\GIS\Selectie_watergangen_21062023\norm_profiles.shp"
output_path = r"D:\Work\Project\HL-23006\GIS\Selectie_watergangen_21062023\sanitized_v3.shp"

geometry_accuracy = 0

branches_gdf = gpd.read_file(branches_path)
duikers_gdf = gpd.read_file(duikers_path).to_crs(branches_gdf.crs)
branches_gdf["STATUSOBJE"] = branches_gdf["STATUSOBJE"].astype(int)
duikers_gdf["STATUSOBJE"] = duikers_gdf["STATUSOBJE"].astype(int)

branches_gdf = branches_gdf.dropna(axis=0, subset="STATUSOBJE")
branches_gdf = branches_gdf.loc[branches_gdf["STATUSOBJE"] == 300, :]
branches_gdf = snap_endpoints(branches_gdf, geometry_accuracy=geometry_accuracy)
duikers_gdf = duikers_gdf.loc[duikers_gdf["STATUSOBJE"] == 300, :]
duikers_gdf = snap_endpoints(duikers_gdf, geometry_accuracy=geometry_accuracy)

duikers_gdf.geometry = duikers_gdf.geometry.buffer(0.5, cap_style=2)
duikers_gdf = gpd.GeoDataFrame(geometry=duikers_gdf.geometry, crs=duikers_gdf.crs)
duikers_gdf = duikers_gdf.dissolve(by=None).explode(ignore_index=True)
# duikers_gdf = gpd.overlay(duikers_gdf, duikers_gdf, how="union")


branches_sel_gdf = gpd.overlay(branches_gdf, duikers_gdf, how="difference").explode(
    ignore_index=True
)
branches_cul_gdf = gpd.overlay(branches_gdf, duikers_gdf, how="intersection").explode(
    ignore_index=True
)
# branches_gdf = snap_endpoints(branches_gdf, geometry_accuracy=geometry_accuracy)


print(branches_sel_gdf.duplicated(subset="CODE").sum())

intersected_gdf = (
    validate_network_topology(branches_sel_gdf, geometry_accuracy=geometry_accuracy)
    .drop(columns=["index_right"])
    .dropna(axis=0, subset="CODE")
)
print(branches_gdf.shape)
print(intersected_gdf.shape)
print(branches_cul_gdf.shape)


comb_gdf = gpd.GeoDataFrame(
    pd.concat([intersected_gdf, branches_cul_gdf]), geometry="geometry", crs=intersected_gdf.crs
)
comb_gdf["globalid"] = [str(uuid.uuid4()) for _ in range(comb_gdf.shape[0])]
comb_gdf.set_index("globalid", inplace=True)
# comb_gdf = snap_nodes(comb_gdf, geometry_accuracy=geometry_accuracy)
comb_gdf.to_file(np_path)
comb_gdf = combine_branches(in_gdf=comb_gdf, geometry_accuracy=geometry_accuracy)
comb_gdf = snap_endpoints(comb_gdf, geometry_accuracy=0)
comb_gdf = combine_branches(in_gdf=comb_gdf, geometry_accuracy=geometry_accuracy)

try:
    os.remove(output_path)
except OSError:
    pass

comb_gdf.to_file(output_path)

print(comb_gdf.shape)

comb_gdf.plot()
plt.show()
