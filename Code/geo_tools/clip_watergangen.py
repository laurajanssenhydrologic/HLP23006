# %% Init
import uuid
from copy import copy
from typing import List, Tuple

import geopandas as gpd
import numpy as np
from scipy.spatial import KDTree
from shapely.geometry import LineString, MultiPoint, Point
from tqdm import tqdm


def check_branch_overlap(
    new_branches: gpd.GeoDataFrame, old_branches: gpd.GeoDataFrame, min_overlap: float = 0.5
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

    # setting index to new_id for easier indexing later
    old_branches.set_index("new_id", inplace=True)
    new_branches.set_index("new_id", inplace=True)

    # find branches that are MultiLineString due to partial intersection
    MLS_branches_bool = new_branches.geometry.type == "MultiLineString"
    # MLS_branches = new_branches.loc[MLS_branches_bool]
    # other_branches = new_branches.loc[~MLS_branches_bool]

    n_split_branches_init = np.sum(MLS_branches_bool)
    # print("Number of split branches: {}".format(n_split_branches_init))

    n_replaced = 0
    n_dropped = 0
    out_branches = copy(new_branches)

    # Loop over all branches. If MLS, check if minimum overlap is achieved and then replace, otherwise, drop.
    # If not MLS, check for sufficient overlap, if not, drop.

    for ix, branch in new_branches.iterrows():
        new_length = np.sum(new_branches.loc[branch.name].geometry.length)
        old_length = np.sum(old_branches.loc[branch.name].geometry.length)

        if MLS_branches_bool[ix]:
            if (new_length / old_length) >= min_overlap:
                out_branches.loc[branch.name] = old_branches.loc[branch.name]
                n_replaced += 1
            else:
                out_branches.drop(index=branch.name, inplace=True)
                n_dropped += 1

        elif not MLS_branches_bool[ix]:
            if (new_length / old_length) < min_overlap:
                out_branches.drop(index=branch.name, inplace=True)
                n_dropped += 1

    n_split_branches = np.sum(out_branches.geometry.type == "MultiLineString")
    # print(
    #     "Replaced {} branches with sufficient overlap. Dropped {} branches without sufficient overlap".format(
    #         n_replaced, n_dropped
    #     )
    # )
    # print("Remaining number of split branches: {}".format(n_split_branches))

    if n_split_branches > 0:
        print(out_branches.loc[out_branches.geometry.type == "MultiLineString"])
        raise ValueError("still found split branches...")

    return out_branches, [n_split_branches_init, n_split_branches, n_replaced, n_dropped]


def clip_branches(
    in_branches_path: str,
    overlay_branches_path: str,
    buffer_dist=5,
    min_overlap: float = 0.5,
    max_distance: float = 0.5,
    min_connectivity: int = 2,
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
    in_branches = (
        gpd.read_file(in_branches_path, geometry="geometry")
        .to_crs(overlay_branches.crs)
        .explode(index_parts=False)
    )
    in_branches["new_id"] = [str(uuid.uuid4()) for _ in range(in_branches.shape[0])]

    # compute extend of input branches and clip overlay branches to it
    convex_hull = in_branches.dissolve(by=None).convex_hull
    overlay_branches_clipped = overlay_branches.clip(convex_hull)

    pbar = tqdm(total=6)
    pbar.set_description("Clipping overlay branches")
    pbar.update(1)

    # buffer overlay branches
    pbar.set_description("Buffering overlay branches")
    buffered_branches = gpd.GeoDataFrame(
        geometry=overlay_branches_clipped.buffer(distance=buffer_dist)
    ).dissolve(by=None)

    pbar.update(1)

    # clip branches using utility function
    pbar.set_description("Clipping input branches with overlay branches")
    out_branches = _clip_branches(in_branches=in_branches, buffered_branches=buffered_branches)

    pbar.update(1)

    # skip branches with an overlap of less than min_overlap
    pbar.set_description("Droping branches with insufficient overlap")
    out_branches, stats = check_branch_overlap(
        new_branches=out_branches, old_branches=in_branches, min_overlap=min_overlap
    )

    pbar.update(1)

    # Validate network toplogy
    pbar.set_description("Validating network toplogy")
    out_branches = validate_network_topology(branches_gdf=out_branches, snap_distance=buffer_dist)

    pbar.update(1)

    # skip branches that are not connected to at least two other branches
    pbar.set_description("Droping branches with too little connectivity")
    out_branches = skip_branches_con(
        in_branches=out_branches, max_distance=max_distance, min_connectivity=min_connectivity
    )

    pbar.update(1)

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


def delete_branches(gdf: gpd.GeoDataFrame, snap_distance: float) -> gpd.GeoDataFrame:
    _gdf = gdf.loc[gdf.geometry.length > snap_distance, :]

    short_branches_ix = np.where(gdf.geometry.length <= snap_distance)[0]
    print(gdf.shape[0], _gdf.shape[0], short_branches_ix)
    boundaries = []
    for branch_ix in short_branches_ix:
        branch = gdf.iloc[branch_ix, :]
        boundary = branch.geometry.boundary
        if boundary not in boundaries:
            boundaries.append(boundary)

    l_points = []
    r_points = []
    for ix, boundary in enumerate(boundaries):
        l_point = Point(boundary.geoms[0].x, boundary.geoms[0].y)
        r_point = Point(boundary.geoms[1].x, boundary.geoms[1].y)

        l_points.append(l_point)
        r_points.append(r_point)

    l_mpoints = MultiPoint(l_points)
    # r_mpoints = MultiPoint(r_points)

    # # print(*boundaries)
    gdf = copy(_gdf)
    for ix, (name, branch) in enumerate(_gdf.iterrows()):
        # Delete z-dimensison of linestring if it exists
        if np.array(branch.geometry.coords).shape[1] > 2:
            geometry = [Point(x, y) for x, y, _ in branch.geometry.coords]
            try:
                gdf.geometry.iloc[ix] = LineString(geometry)
            except ValueError:
                print("oepsie")
        else:
            geometry = [Point(x, y) for x, y in branch.geometry.coords]

        boundary = branch.geometry.boundary

        if l_mpoints.contains(geometry[0]):
            jx = l_points.index(geometry[0])
            _branch = copy(branch)
            geometry[0] = r_points[jx]
            _branch.geometry = LineString(geometry)
            gdf.iloc[ix] = _branch

        if l_mpoints.contains(geometry[-1]):
            jx = l_points.index(geometry[-1])
            _branch = copy(branch)
            geometry[-1] = r_points[jx]
            _branch.geometry = LineString(geometry)
            gdf.iloc[ix] = _branch

    return gdf


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


def skip_branches_con(
    in_branches: gpd.GeoDataFrame, max_distance: float = 0.5, min_connectivity: int = 2
) -> gpd.GeoDataFrame:
    """
    Function that removes branches from GeoDataFrame if they are not sufficiently connected to other branches
    The criterion for this is min_connectivity for the number of branches that should coincide in a point
    Max_distance is the distance tolerance for branches starting and ending in a point

    Args:
        in_branches (gpd.GeoDataFrame): GeoDataFrame with branches containging a geometry column with LineStrings
        max_distance (float): maximum distance at which branches are considered connected
        min_connectivity (int): minumum number of branches that should meet in a point

    Returns:
        out_branhces (gpd.GeoDataFrame): GeoDataFrame with branches that fullfill the connectivity criterion
    """

    # print("Checking for isolated branches")
    n_in = in_branches.shape[0]

    # skip branches further than max_dist away from another branch
    point_list = []
    for ix, branch in in_branches.iterrows():
        coords = branch.geometry.coords
        point_list.append(coords[0])
        point_list.append(coords[-1])

    # build spatial KDTree from start and end points
    kdtree = KDTree(point_list)
    sdm = kdtree.sparse_distance_matrix(
        kdtree, max_distance=max_distance, output_type="coo_matrix"
    )

    # matrix is symetrical, so picking one of the two axis to count number of non empty entries (explicit zeros allowed)
    nnz = sdm.getnnz(axis=0)

    # compare number of points within max_distance to min_connectivity
    bool_array = nnz < min_connectivity

    # Drop branches with insufficient connectivity on both sides
    out_branches = copy(in_branches)
    for ix, (name, branch) in enumerate(in_branches.iterrows()):
        if (bool_array[ix * 2]) & (bool_array[ix * 2 + 1]):
            out_branches.drop(index=name, inplace=True)
        # elif bool_array[ix * 2 + 1]:
        #     out_branches.drop(index=name, inplace=True)

    # out_branches = in_branches
    n_out = out_branches.shape[0]
    # print("Dropped {} isolated branches".format(n_in - n_out))
    return out_branches


def validate_network_topology(
    branches_gdf: gpd.GeoDataFrame, snap_distance: float = None
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
    buffered_branches.geometry = buffered_branches.geometry.buffer(0.1)
    gdf = gpd.GeoDataFrame(union_result, columns=["geometry"], geometry="geometry", crs=28992)

    # add data from buffered branches if lines in gdf fall within. But keep geometry of branches in gdf
    intersected_gdf = gdf.sjoin(buffered_branches, how="left", predicate="within")

    if snap_distance is not None:
        if snap_distance > 0:
            intersected_gdf = delete_branches(gdf=intersected_gdf, snap_distance=snap_distance)
        else:
            raise ValueError("snap_distance should be > 0")

    return intersected_gdf


# %% Initialize
print("initialize")
p_folder = r"D:\work\P1414_ROI\GIS"
version = "v10"
# p_folder = r"D:\work\Project\P1414\GIS"
old_rm_branches_path = p_folder + r"\Randstadmodel_oud\rm_Branches_28992_edited_v12.shp"

agv_branches_path = p_folder + r"\WAGV\hydroobject_v13\hydroobject_v13_clipped.shp"
clipped_agv_branches_path = p_folder + r"\Uitgesneden watergangen\AGV_{}_test.shp".format(version)

HDSR_branches_path = p_folder + r"\HDSR\Legger\Hydro_Objecten(2)\HydroObject.shp"
clipped_HDSR_branches_path = p_folder + r"\Uitgesneden watergangen\HDSR_{}_test.shp".format(
    version
)

HHD_branches_path = (
    p_folder + r"\HHDelfland\Legger_Delfland_shp\Oppervlaktewaterlichamen\Primair water.shp"
)
clipped_HHD_branches_path = p_folder + r"\Uitgesneden watergangen\HHD_{}_test.shp".format(version)

HHR_branches_path = p_folder + r"\HHRijnland\Niet legger\Watergangen_edited_v3.shp"
clipped_HHR_branches_path = p_folder + r"\Uitgesneden watergangen\HHR_{}_test.shp".format(version)

HHR_west_branches_path = p_folder + r"\HHRijnland\Niet legger\Westeinderplassen_cut_v5.shp"
clipped_HHR_west_branches_path = (
    p_folder + r"\Uitgesneden watergangen\HHR_west_{}_test.shp".format(version)
)

HHSK_branches_path = p_folder + r"\HHSK\Legger\Hoofdwatergang.shp"
clipped_HHSK_branches_path = p_folder + r"\Uitgesneden watergangen\HHSK_{}_test.shp".format(
    version
)
RMM_branches_path = p_folder + r"\Rijn Maasmonding\RMM_Branches_edited.shp"
fixed_RMM_branches_path = p_folder + r"\Uitgesneden watergangen\RMM_branches.shp"

AGV = False
HDSR = False
HHD = False
HHR = False
HHSK = False
RMM = False
HHR_westeinder = True

# %% AGV
if AGV:
    print("AGV")

    intersected_branches = clip_branches(
        in_branches_path=agv_branches_path,
        overlay_branches_path=old_rm_branches_path,
    )
    intersected_branches.to_file(clipped_agv_branches_path)

# %% HDSR
if HDSR:
    print("HDSR")

    intersected_branches = clip_branches(
        in_branches_path=HDSR_branches_path,
        overlay_branches_path=old_rm_branches_path,
    )
    intersected_branches.to_file(clipped_HDSR_branches_path)

# %% HHD
if HHD:
    print("HHD")

    intersected_branches = clip_branches(
        in_branches_path=HHD_branches_path, overlay_branches_path=old_rm_branches_path
    )
    intersected_branches.to_file(clipped_HHD_branches_path)

# %% HHR
if HHR:
    print("HHR")

    intersected_branches = clip_branches(
        in_branches_path=HHR_branches_path, overlay_branches_path=old_rm_branches_path
    )
    intersected_branches.to_file(clipped_HHR_branches_path)

# %% HHSK
if HHSK:
    print("HHSK")

    intersected_branches = clip_branches(
        in_branches_path=HHSK_branches_path,
        overlay_branches_path=old_rm_branches_path,
    )
    intersected_branches.to_file(clipped_HHSK_branches_path)

if RMM:
    print("RMM")
    in_branches = gpd.read_file(RMM_branches_path, geometry="geometry", crs=28992).explode(
        index_parts=False
    )
    in_branches["new_id"] = [str(uuid.uuid4()) for _ in range(in_branches.shape[0])]
    out_branches = validate_network_topology(branches_gdf=in_branches, snap_distance=0.5)
    out_branches.to_file(fixed_RMM_branches_path)

# %% HHR Westeinderplassen
if HHR_westeinder:
    print("HHR Westeinder")

    intersected_branches = clip_branches(
        in_branches_path=HHR_west_branches_path, overlay_branches_path=old_rm_branches_path
    )
    intersected_branches.to_file(clipped_HHR_west_branches_path)
