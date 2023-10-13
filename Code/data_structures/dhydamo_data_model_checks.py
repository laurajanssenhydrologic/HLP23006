import uuid

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import LineString, Point


def geometry_check(geometry: gpd.GeoSeries) -> bool:
    """
    Function that checks if provided geometry is a shapely Point or LineString

    Args:
        geometry(gpd.GeoSeries): input geometry to check

    Returns:
        results (bool): True, if geometry is a Point or LineString; False if not
    """
    # Allow points and linestrings
    if isinstance(geometry, Point) or isinstance(geometry, LineString):
        return True
    # # Verify that multipoints are infact normal points
    # elif isinstance(geometry, MultiPoint):
    #     if len(geometry.geoms) == 1:
    #         return True
    #     else:
    #         return False
    else:
        return False


def globalid_check(globalid: str) -> bool:
    """
    Function that checks if provided globaid is a UUID

    Args:
        globalid (str): globalid to check

    Returns:
        results (bool): True, if globalid is a UUID; False if not
    """
    try:
        uuid_obj = uuid.UUID(globalid)
        return True
    except ValueError:
        return False


def none_geometry_check(geometry: gpd.GeoSeries) -> bool:
    """
    Function that checks if provided geometry is None

    Args:
        geometry(gpd.GeoSeries): input geometry to check

    Returns:
        results (bool): True, if geometry is None; False if not
    """
    # Allow only None types
    if geometry is None:
        return True
    else:
        return False


def validate_codes(values: dict, struct_list: list = ["brug", "duiker", "gemaal", "stuw"]) -> dict:

    codes = None
    # check if field is assigned and if so add code column to list of codes
    for ix, struct in enumerate(struct_list):
        if values.get(struct) is not None:
            if codes is None:
                codes = values.get(struct)["code"]
            else:
                codes = pd.concat([codes, values.get(struct)["code"]], ignore_index=True)

    # if code columns were assigned, check for duplicates
    if codes is not None:
        duplicate_codes = codes.duplicated(keep=False)
        if np.sum(duplicate_codes) > 0:
            print("The codes that cause the errors are:")
            print(codes.loc[duplicate_codes])
            raise ValueError("Duplicate codes found in structures")
    return values
