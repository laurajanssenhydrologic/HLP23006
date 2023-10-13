# %% Init
import uuid
from copy import copy
from typing import List, Tuple, Union

import geopandas as gpd
import networkx as nx
import numpy as np
from scipy.spatial import KDTree
from shapely.geometry import LineString, MultiPoint, Point
from tqdm import tqdm

from network_tools import combine_straight_branches
from networkx_tools import gdf_to_nx


def check_branch_overlap(
    new_branches: gpd.GeoDataFrame,
    old_branches: gpd.GeoDataFrame,
    geometry_accuracy: float,
    min_overlap: float,
    max_graphs: int,
) -> Tuple[gpd.GeoDataFrame, List[float]]:
    """
    Computes overlap between unclipped and clipped branches.
    Only keeps branches with an overlap larger than min_overlap.
    Replaces MultiLineStrings (i.e. lines that are only partiallyt clipped) with unclipped data if overlap is larger than min_ocerlap

    Args:
        new_branches (gpd.GeoDataFrame): clipped branches that are checked against old_branches
        old_branches (gpd.GeoDataFrame): un-clipped branches, used to compute the amount of overlap between clipped and unclipped
        min_overlap (float): fraction of clipped branch / unclipped branch. Branches are kept if this fraction is larger than min_overlap

    Returns:
        out_branches (gpd.GeoDataFrame): sanitized branches
        stats (List): list of stats [number of initial split branches, number of remaining split branches,
                                        number of replaced branches, number of dropped branches]

    """

    # print("Checking for split branches and fixing them")

    # G = gdf_to_nx(old_branches)
    _out_branches = copy(old_branches)
    _out_branches["new_length"] = 0
    _out_branches["old_length"] = 0
    overlap_list = []
    for ix, branch in new_branches.iterrows():
        new_length = np.sum(
            new_branches[new_branches["new_id"] == branch["new_id"]].geometry.length
        )
        old_length = np.sum(
            old_branches[old_branches["new_id"] == branch["new_id"]].geometry.length
        )
        overlap = new_length / old_length

        if np.isnan(new_length):
            new_length = 0

        _out_branches.loc[_out_branches["new_id"] == branch["new_id"], "new_length"] = new_length
        _out_branches.loc[_out_branches["new_id"] == branch["new_id"], "old_length"] = old_length

        overlap_list.append({"id": ix, "overlap": overlap})

    G = gdf_to_nx(_out_branches, decimals=geometry_accuracy)
    component_list = sorted(nx.connected_components(G), key=len, reverse=True)
    H = None
    for n in range(min(max_graphs, len(component_list))):
        if n > len(component_list):
            break

        S = G.subgraph(component_list[n])

        new_lengths = old_lengths = 0
        for x, y, data in S.edges(data=True):
            if data["new_length"] > 0:
                new_lengths += data["new_length"]
                old_lengths += data["geometry"].length
        if new_lengths > 1e3:
            overlap = new_lengths / old_lengths
        else:
            overlap = 0
        print(overlap)
        if overlap >= min_overlap:
            if H is None:
                H = nx.Graph(S)
            else:
                H = nx.compose(H, nx.Graph(S))

    H = combine_straight_branches(G=H)

    out_branches_list = []
    for x, y, data in H.edges(data=True):
        out_branches_list.append(data)

    out_branches = gpd.GeoDataFrame(
        data=out_branches_list, geometry="geometry", crs=new_branches.crs
    )

    return out_branches, None


def clip_branches(
    in_branches_path: str,
    overlay_branches_path: str,
    buffer_dist=10,
    geometry_accuracy: int = 0,
    in_branches: gpd.GeoDataFrame = None,
    max_graphs: float = 100,
    min_overlap: float = 0.25,
    min_width: float = 5,
    width_column: str = None,
) -> gpd.GeoDataFrame:
    """
    Wrapper function that performs all actions required to clip branches using an overlay polygon
    and drops branches that are not sufficiently overlaid or not connected to at least two other branches

    Args:
        in_branches_path (str): path where input branches shape is stored
        overlay_branches_path (str): path where overlay shape is stored
        buffer_dist (float): distance around lines where buffer is applied
        min_overlap (float): fraction (0-1) of branch that should be overlaped by buffer
        max_distance (float): max_distance that start and end of branch are allowed to lay apart
                              used to check how many branches are connected in one point
        min_connectivity (int): number indicating how many branches should connect in one point in order to keep a branch
                                default 2 means 1 other branch should connect

    Returns:
        out_branches (gpd.GeoDataFrame): GeoDataFrame contianing the clipped branches that fullfill two criterion
                                         1. fraction of branch that is covered by overlay branches is >= min_overlap
                                         2. at least 1 other branch connects to both ends of a branch

    """

    # Load input branches, explode multilinestrings, and assign unique ids.
    # Also load overlay branches
    overlay_branches, _ = read_rm_branches(overlay_branches_path)
    if in_branches is None:
        if in_branches_path.endswith(r".shp"):
            in_branches = (
                gpd.read_file(in_branches_path, geometry="geometry")
                .to_crs(overlay_branches.crs)
                .explode(index_parts=False)
            )
        elif in_branches_path.endswith(r".gpkg"):
            in_branches = (
                gpd.read_file(in_branches_path, layer="waterloop", geometry="geometry")
                .to_crs(overlay_branches.crs)
                .explode(index_parts=False)
            )
    in_branches["new_id"] = [str(uuid.uuid4()) for _ in range(in_branches.shape[0])]
    if (width_column is not None) and (min_width is not None):
        in_branches = in_branches.loc[in_branches[width_column].astype(float) > min_width, :]

    # compute extend of input branches and clip overlay branches to it
    convex_hull = in_branches.dissolve(by=None).convex_hull
    overlay_branches_clipped = overlay_branches.clip(convex_hull)

    pbar = tqdm(total=5)
    pbar.set_description("Clipping overlay branches")
    pbar.update(1)

    # buffer overlay branches
    pbar.set_description("Buffering overlay branches")
    buffered_branches = gpd.GeoDataFrame(
        geometry=overlay_branches_clipped.buffer(distance=buffer_dist)
    ).dissolve(by=None)

    pbar.update(1)

    # Validate network toplogy
    pbar.set_description("Validating network toplogy")
    out_branches = validate_network_topology(
        branches_gdf=in_branches, geometry_accuracy=geometry_accuracy
    )
    in_branches = copy(out_branches)
    in_branches["new_id"] = [str(uuid.uuid4()) for _ in range(in_branches.shape[0])]

    pbar.update(1)

    # clip branches using utility function
    pbar.set_description("Clipping input branches with overlay branches")
    out_branches = _clip_branches(in_branches=in_branches, buffered_branches=buffered_branches)

    pbar.update(1)

    # skip branches with an overlap of less than min_overlap
    pbar.set_description("Droping branches with insufficient overlap")
    out_branches, stats = check_branch_overlap(
        geometry_accuracy=geometry_accuracy,
        max_graphs=max_graphs,
        min_overlap=min_overlap,
        new_branches=out_branches,
        old_branches=in_branches,
    )

    # # skip branches that are not connected to at least two other branches
    # pbar.set_description("Droping branches with too little connectivity")
    # out_branches = skip_branches_con(
    #     in_branches=out_branches, max_distance=max_distance, min_connectivity=min_connectivity
    # )

    # pbar.update(1)

    return out_branches


def _clip_branches(
    in_branches: gpd.GeoDataFrame,
    buffered_branches: gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
    """
    Utility function that clippes the linestrings in in_branches using polygon buffered_branches.

    Args:
        in_branches (gpd.GeoDataFrame): input GeoDataFrame with branches as linestrings
        buffered_branches (gpd.GeoDataFrame): GeoDataFrame containing one polygon to use as overlay mask

    Returns:
        out_branches (gpd.GeoDataFrame): output GeoDataFrame with branches as linestrings
    """
    if buffered_branches.shape[0] > 1:
        raise ValueError(
            "Clipping with more than one polygon might result in duplications in the results\nplease use only one polygon for clipping"
        )

    # print("Intersecting buffered branches with input branches")
    # intersect branches with buffered branches from old model
    out_branches = gpd.overlay(
        in_branches, buffered_branches, how="intersection", keep_geom_type=True
    )

    # return out_branches
    return out_branches


def read_rm_branches(
    rm_branches_path: str, epsg: int = 28992
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """
    function to read old RM-model branches and skip underpasses (marked by "OND")

    Args:
        rm_branches_path (str): path of shape file containing the old branches

    Returns:
        rm_branches (gpd.GeoDataFrame): GeoDataFrame containing the branches
        onderdoorgangen (gpd.GeoDataFrame): GeoDataFrame containing underpasses
    """

    rm_branches = gpd.read_file(rm_branches_path, geometry="geometry").to_crs(crs=epsg)
    ondd_bool = (rm_branches["Source"].str.contains("OND")) | (
        rm_branches["Target"].str.contains("OND")
    )
    onderdoorgangen = rm_branches.loc[ondd_bool, :]
    rm_branches = rm_branches.loc[~ondd_bool, :]
    return rm_branches, onderdoorgangen


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


# %% Initialize
if __name__ == '__main__':
    print("initialize")
    p_folder = r"D:\Work\Project\P1414\GIS"
    # p_folder = r"D:\work\P1414_ROI\GIS"
    version = "v3"
    # p_folder = r"D:\work\Project\P1414\GIS"
    old_rm_branches_path = p_folder + r"\Randstadmodel_oud\rm_Branches_28992_edited_v10.shp"

    # agv_branches_path = p_folder + r"\WAGV\hydroobject_v13\hydroobject_v13_clipped.shp"
    agv_branches_path = p_folder + r"\WAGV\hydrovak\hydrovak.shp"
    clipped_agv_branches_path = p_folder + r"\Uitgesneden watergangen\AGV_{}.shp".format(version)

    HDSR_branches_path = p_folder + r"\HDSR\Legger\Hydro_Objecten(2)\HydroObject.shp"
    clipped_HDSR_branches_path = p_folder + r"\Uitgesneden watergangen\HDSR_{}.shp".format(version)

    HHD_branches_path = (
        p_folder + r"\HHDelfland\Legger_Delfland_shp\Oppervlaktewaterlichamen\Primair water.shp"
    )
    clipped_HHD_branches_path = p_folder + r"\Uitgesneden watergangen\HHD_{}.shp".format(version)

    HHR_branches_path = p_folder + r"\HHRijnland\Legger\Watergang\Watergang_as.shp"
    clipped_HHR_branches_path = p_folder + r"\Uitgesneden watergangen\HHR_{}.shp".format(version)

    # HHR_west_branches_path = p_folder + r"\HHRijnland\Niet legger\Westeinderplassen_cut_v5.shp"
    # clipped_HHR_west_branches_path = (
    #     p_folder + r"\Uitgesneden watergangen\HHR_west_{}_test.shp".format(version)
    # )

    # HHSK_branches_path = p_folder + r"\HHSK\Legger\Hoofdwatergang.shp"
    HHSK_branches_path = p_folder + r"\HHSK\Legger\oppervlaktewater_lijnen.shp"
    # HHSK_branches_path = p_folder + r"\HHSK\Legger\Hoofdwatergang_nieuw\Hoofdwater_primair.shp"
    clipped_HHSK_branches_path = p_folder + r"\Uitgesneden watergangen\HHSK_{}.shp".format(version)
    # clipped_HHSK_branches_path = p_folder + r"\Uitgesneden watergangen\HHSK_{}_test.shp".format(
    #     version
    # )
    RMM_branches_path = p_folder + r"\Rijn Maasmonding\without_lek\RMM_Branches.shp"
    fixed_RMM_branches_path = p_folder + r"\Rijn Maasmonding\without_lek\RMM_Branches_fixed.shp"

    AGV = False
    HDSR = False
    HHD = True
    HHR = True
    HHSK = True
    RMM = False
    HHR_westeinder = False

    max_graphs = 100
    min_width = 0
    min_overlap = 0.25

    # %% AGV
    if AGV:
        print("AGV")

        intersected_branches = clip_branches(
            in_branches_path=agv_branches_path,
            overlay_branches_path=old_rm_branches_path,
            max_graphs=max_graphs,
            min_width=min_width,
            min_overlap=min_overlap,
            width_column="IWS_W_WATB",
        )
        intersected_branches.to_file(clipped_agv_branches_path)

    # %% HDSR
    if HDSR:
        print("HDSR")
        branches = (
            gpd.read_file(HDSR_branches_path, geometry="geometry")
            .to_crs(28992)
            .explode(index_parts=False)
        )
        branches = branches.loc[branches["CATEGORIEO"] == 1, :]
        intersected_branches = clip_branches(
            in_branches_path=HDSR_branches_path,
            overlay_branches_path=old_rm_branches_path,
            geometry_accuracy=0,
            in_branches=branches,
            max_graphs=max_graphs,
            min_width=min_width,
            min_overlap=min_overlap,
            width_column="IWS_W_WATB",
        )
        intersected_branches.to_file(clipped_HDSR_branches_path)

    # %% HHD
    if HHD:
        print("HHD")
        branches = (
            gpd.read_file(HHD_branches_path, geometry="geometry")
            .to_crs(28992)
            .explode(index_parts=False)
        )
        branches = branches.loc[branches["FUNCTIE"] == "Primair boezemwater", :]
        intersected_branches = clip_branches(
            in_branches_path=HHD_branches_path,
            overlay_branches_path=old_rm_branches_path,
            in_branches=branches,
            max_graphs=max_graphs,
            min_width=None,
            min_overlap=min_overlap,
            width_column=None,
        )
        intersected_branches.to_file(clipped_HHD_branches_path)

    # %% HHR
    if HHR:
        print("HHR")
        branches = (
            gpd.read_file(HHR_branches_path, geometry="geometry")
            .to_crs(28992)
            .explode(index_parts=False)
        )
        branches = branches.loc[branches["CATEGORIEO"] == "primair", :]
        intersected_branches = clip_branches(
            in_branches_path=HHR_branches_path,
            overlay_branches_path=old_rm_branches_path,
            in_branches=branches,
            max_graphs=max_graphs,
            min_width=min_width,
            min_overlap=min_overlap,
            width_column="BREEDTE",
        )
        intersected_branches.to_file(clipped_HHR_branches_path)

    # %% HHSK
    if HHSK:
        print("HHSK")
        branches = (
            gpd.read_file(HHSK_branches_path, geometry="geometry")
            .to_crs(28992)
            .explode(index_parts=False)
        )
        branches = branches.loc[branches["STATUSOBJE"] == 3, :]
        branches = branches.loc[branches["CATEGORIEO"] == 1, :]
        intersected_branches = clip_branches(
            in_branches_path=HHSK_branches_path,
            overlay_branches_path=old_rm_branches_path,
            in_branches=branches,
            max_graphs=max_graphs,
            min_width=min_width,
            min_overlap=min_overlap,
            width_column="BREEDTE",
        )
        intersected_branches.to_file(clipped_HHSK_branches_path)

    if RMM:
        print("RMM")
        in_branches = gpd.read_file(RMM_branches_path, geometry="geometry", crs=28992).explode(
            index_parts=False
        )
        in_branches["new_id"] = [str(uuid.uuid4()) for _ in range(in_branches.shape[0])]
        out_branches = validate_network_topology(branches_gdf=in_branches, geometry_accuracy=0)
        out_branches.to_file(fixed_RMM_branches_path)

# %%
