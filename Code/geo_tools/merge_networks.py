import numpy as np
import shapely
from scipy.spatial import KDTree


def non_zero_min_idx(sdm_array):
    """
    This function returns the index with the smallest distance from an array taken from the
    sparse distance matrix
    """
    non_zeros = np.nonzero(sdm_array)[0]
    try:
        min_idx = non_zeros[0]
        for i in non_zeros:
            if sdm_array[i] < sdm_array[min_idx]:
                min_idx = i

    except:
        min_idx = None

    return min_idx, sdm_array[min_idx]


def check_more_branches(sdm_array, base_idx, min_idx):
    """
    This function checks whether there are more branches connecting to a single branch, and
    returns whether the given index is the closest branch
    """
    # TODO: Repair when needed, not functioning properly

    if np.sum(sdm_array[:, min_idx]) == 1:
        check_min_idx = True

    else:
        min_match_idx, _ = non_zero_min_idx(sdm_array[:, min_idx])

        if min_match_idx == base_idx:
            check_min_idx = True
        else:
            check_min_idx = False

    return check_min_idx


def merge_networks(data_base_input, data_match_input, max_dist=1, outputfile_path=None):
    """
    This function snaps geometries from the base dataset to the matching dataset if they are within
    a specified distance from eachother.

    Arguments:
        - data_base_input (GeoDataFrame): The base geodataframe that gets snapped
        - data_match_input (GeoDataFrame): The geodataframe to which the geometries from the base data will be snapped
        - max_dist (int): The maximum distance within which branches are snapped to eachother
        - outputfile_path (str): Location where to save the resulting dataset

    Returns:
        - data_base (GeoDataFrame): A copy of the original data_base_input, but with the geometries
                                    adjusted to the data_match_input
    """
    data_base = data_base_input.copy()
    data_match = data_match_input.copy()

    point_list_base = []
    for ix, branch in data_base.iterrows():
        coords = branch.geometry.coords
        x1, y1, *_ = coords[0]
        x2, y2, *_ = coords[-1]
        point_list_base.append((x1, y1))
        point_list_base.append((x2, y2))

    point_list_match = []
    for ix, branch in data_match.iterrows():
        coords = branch.geometry.coords
        x1, y1, *_ = coords[0]
        x2, y2, *_ = coords[-1]
        point_list_match.append((x1, y1))
        point_list_match.append((x2, y2))

    # build spatial KDTree from start and end points
    kdtree_base = KDTree(point_list_base)
    kdtree_match = KDTree(point_list_match)

    # Set up a sparse distance matrix
    sdm = kdtree_base.sparse_distance_matrix(
        kdtree_match, max_distance=max_dist, output_type="coo_matrix"
    )

    # Check for non-zero (explicit zeroes allowed) in base dataset
    nnz_base = sdm.getnnz(axis=1)
    bool_array = nnz_base > 0

    # Get the geometries of the branches
    geometry_branches = data_base.geometry

    for ix, (name, branch) in enumerate(data_base.iterrows()):
        if bool_array[ix * 2] and bool_array[ix * 2 + 1]:

            min_idx_start, min_dist_start = non_zero_min_idx(sdm.toarray()[ix * 2])
            min_idx_end, min_dist_end = non_zero_min_idx(sdm.toarray()[ix * 2 + 1])

            if min_dist_start <= min_dist_end:
                min_idx_coords = min_idx_start
                coord_index = ix * 2
                replace_coord = 0
            else:
                min_idx_coords = min_idx_end
                coord_index = ix * 2 + 1
                replace_coord = -1
            try:
                closest_branch_check = check_more_branches(
                    sdm.toarray(), coord_index, min_idx_coords
                )
            except:
                min_idx_coords = None

            if min_idx_coords != None and closest_branch_check:
                coords = branch.geometry.coords[:]
                # print(min_idx_coords)
                coords[replace_coord] = point_list_match[min_idx_coords]
                geometry_branches[ix] = shapely.geometry.LineString(coords)

        elif bool_array[ix * 2]:

            min_idx_coords, _ = non_zero_min_idx(sdm.toarray()[ix * 2])

            if min_idx_coords != None:
                coords = branch.geometry.coords[:]
                coords[0] = point_list_match[min_idx_coords]
                geometry_branches[ix] = shapely.geometry.LineString(coords)

        elif bool_array[ix * 2 + 1]:

            min_idx_coords, _ = non_zero_min_idx(sdm.toarray()[ix * 2 + 1])

            if min_idx_coords != None:
                coords = branch.geometry.coords[:]
                coords[-1] = point_list_match[min_idx_coords]
                geometry_branches[ix] = shapely.geometry.LineString(coords)

    data_base = data_base.set_geometry(geometry_branches)

    if outputfile_path is not None:
        data_base.to_file(outputfile_path)

    return data_base
