import uuid
from copy import copy

import geopandas as gpd
import networkx as nx
import numpy as np
from shapely.geometry import LineString, MultiLineString
from shapely.ops import linemerge


def check_for_nan(value):
    nan = False
    if value is None:
        nan = True
    elif isinstance(value, (int, float)):
        if np.isnan(value):
            nan = True
    return nan


def combine_straight_branches(G: nx.Graph):
    # Select all nodes with only 2 neighbors
    nodes_to_remove = [n for n in G.nodes if len(list(G.neighbors(n))) == 2]

    for node in nodes_to_remove:
        edges = G.edges(node, data=True)
        if len(edges) != 2:
            continue
        elif len(list(G.neighbors(node))) != 2:
            continue
        else:
            try:
                (_, _, data1), (_, _, data2) = edges
            except Exception as e:
                print(edges)
                print(len(edges))
                # print(G(node))
                # raise e
                print(e)
                continue

        data_new = {}
        succesfull = True
        for key, value in data1.items():
            try:
                if check_for_nan(value) and check_for_nan(data2[key]):
                    data_new[key] = None
                    continue
                elif check_for_nan(value):
                    data_new[key] = data2[key]
                    continue
                elif check_for_nan(data2[key]):
                    data_new[key] = value
                    continue

                if isinstance(value, str):
                    data_new[key] = value

                elif isinstance(value, (int, float)):
                    data_new[key] = np.nanmean([value, data2[key]])

                elif isinstance(value, LineString):
                    line = linemerge([data1[key], data2[key]])
                    if isinstance(line, MultiLineString):
                        succesfull = False
                        break
                    else:
                        data_new[key] = line

                elif value is None:
                    data_new[key] = None

                else:
                    print(data1["CODE"])
                    print(value)
                    print(type(value))
                    raise TypeError("Unimplemented type")
            except Exception as e:
                print(data1[key])
                print(data2[key])
                print(e)
                data_new[key] = None

        if succesfull:
            G.add_edge(*G.neighbors(node), **data_new)
            G.remove_node(node)
    return G


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

    def snap_nodes(in_branches: gpd.GeoDataFrame, geometry_accuracy: float):
        in_branches["t_id"] = [str(uuid.uuid4()) for _ in range(in_branches.shape[0])]
        in_branches.set_index("t_id", inplace=True)
        out_branches = copy(in_branches)
        for ix, branch in in_branches.iterrows():

            point_list = []
            for x, y, *_ in branch.geometry.coords[:]:
                point_list.append(
                    (np.around(x, geometry_accuracy), np.around(y, geometry_accuracy))
                )

            geometry = LineString(point_list)
            if geometry.length > 0:
                out_branches.loc[ix, "geometry"] = LineString(point_list)
            else:
                out_branches = out_branches.drop(index=ix)

        return out_branches

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
    buffered_branches.geometry = buffered_branches.geometry.buffer(0.1)
    gdf = gpd.GeoDataFrame(union_result, columns=["geometry"], geometry="geometry", crs=28992)

    # add data from buffered branches if lines in gdf fall within. But keep geometry of branches in gdf
    intersected_gdf = gdf.sjoin(buffered_branches, how="left", predicate="within")

    intersected_gdf = snap_nodes(in_branches=intersected_gdf, geometry_accuracy=geometry_accuracy)
    # if snap_distance is not None:
    #     if snap_distance > 0:
    #         intersected_gdf = delete_branches(gdf=intersected_gdf, snap_distance=snap_distance)
    #     else:
    #         raise ValueError("snap_distance should be > 0")

    return intersected_gdf
