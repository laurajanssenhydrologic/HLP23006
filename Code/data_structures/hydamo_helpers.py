import importlib
import uuid
import warnings
from copy import copy
from typing import Any, List, Tuple

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely import affinity
from shapely.geometry import LineString, MultiPoint, Point

from data_structures.config_templates import FMGisFile
from data_structures.hydamo_globals import (
    BRANCH_FRICTION_FUNCTION,
    MANAGEMENT_DEVICE_TYPES,
    ROUGHNESS_MAPPING_LIST,
    WEIR_MAPPING,
)
from data_structures.roi_data_model import FMDataModel as Datamodel

warnings.filterwarnings(action="ignore", message="Mean of empty slice")


def check_and_fix_duplicate_code(gdf: gpd.GeoDataFrame, column="code") -> gpd.GeoDataFrame:
    """
    Function that takes a GeoDataFrame and a column name (default: column="code"),
    checks if the column exists, and if it does, it checks if there are any duplicate values.
    If there are, it pads the duplicate values with a number

    Args:
      gdf (gpd.GeoDataFrame): gpd.GeoDataFrame
      column (str): the name of the column to check for duplicates. Defaults to code

    Returns:
      gdf(gpd.GeoDataFrame): A GeoDataFrame with the duplicated codes fixed.
    """

    if column in gdf.columns:
        gdf[column] = gdf[column].astype("str").str.strip()
        # _gdf = copy(gdf)
        duplicates = gdf.duplicated(subset=column, keep=False)
        duplicate_codes = gdf.loc[duplicates, column]

        for code in duplicate_codes.unique():
            n_duplicates = np.sum(gdf[column] == code)
            pad_list = [r"{}_{}".format(code, n) for n in np.arange(n_duplicates, dtype=np.int8)]

            gdf.loc[gdf[column] == code, column] = pad_list

    return gdf


def check_column_is_numerical(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Helper function that tries to convert the column in the geodataframe to int or float

    Args:
      gdf (gpd.GeoDataFrame): the GeoDataFrame that you want to check

    Returns:
      gdf (gpd.GeoDataFrame): the checked GeoDataFrame
    """

    if (not isinstance(gdf, int)) | (not isinstance(gdf, float)):
        gdf = gdf.astype(float).fillna(value=np.nan)
    return gdf


def check_is_not_na_number(input: Any, zero_allowed=False) -> bool:
    """
    Helper funciton that checks if the input is not None and not NaN and (if zero_allowed is False) not zero.
    Returns True if this is the case, and false if any of these values are found.

    Args:
      input (Any): The input value to check.
      zero_allowed (bool): If the input is zero, should it be considered a valid number?. Defaults to False

    Returns:
      result (bool): whether the number is na, or valid
    """
    if (input is None) or (np.isnan(input)):
        return False

    if (not zero_allowed) and (input == 0):
        return False

    return True


def check_roughness(structure: gpd.GeoSeries, rougness_map: List = ROUGHNESS_MAPPING_LIST) -> str:
    """
    Function that takes a GeoSeries of structures and a list of roughness types,
    and returns the roughness type  of the structure as a str

    Args:
      structure (gpd.GeoSeries): gpd.GeoSeries
      rougness_map (List): Roughness types with indices

    Returns:
      type_ruwheid (str): the roughenss type as a strinng
    """
    type_ruwheid = structure["typeruwheid"]
    if isinstance(type_ruwheid, int) or isinstance(type_ruwheid, float):
        type_ruwheid = rougness_map[int(type_ruwheid) - 1]
    return type_ruwheid

def create_streefpeil(gdf_peil: gpd.GeoDataFrame, gdf_gemaal: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Function that overlays the peilgebieden shape with the gemalen shape  to determine the streefpeil for each pump
    """
    gdf_peil_crs=gdf_peil.to_crs(gdf_gemaal.crs)
    overlay_peil=gpd.sjoin(gdf_gemaal,gdf_peil_crs,how='inner',predicate='intersects')
    gdf_gemaal_peil=overlay_peil[['code','globalid','geometry','maximalecapaciteit','doelvariabele','peil_marge','ZOMERPEIL','WINTERPEIL','BOVENPEIL','ONDERPEIL','VASTPEIL','streefwaarde']]

    for ix, pump in gdf_gemaal_peil.iterrows():
        #als er een vastpeil is, dan nemen we die
        #anders nemen we het bovenpeil
        #en anders nemen we zomerpeil voor nu
        if np.isnan(pump["VASTPEIL"]) and ~np.isnan(pump["BOVENPEIL"]):
            gdf_gemaal_peil.loc[ix,'streefwaarde']=pump["BOVENPEIL"] 
        elif np.isnan(pump["VASTPEIL"]) and np.isnan(pump["BOVENPEIL"]):
            gdf_gemaal_peil.loc[ix,'streefwaarde']=pump["ZOMERPEIL"]
        else: 
            gdf_gemaal_peil.loc[ix,'streefwaarde']=pump["VASTPEIL"]

    #keep rows with lowest streefwaarde
    gdf_gemaal_peil=gdf_gemaal_peil.sort_values(by='streefwaarde')      
    gdf_gemaal_v2=gdf_gemaal_peil.drop_duplicates(subset='globalid',keep='first')
    gdf_gemaal_v2=gdf_gemaal_v2.sort_index()
    return gdf_gemaal_v2

def create_brugprofiel(
        branches_gdf: gpd.GeoDataFrame,
        dist_tol: float = 0.25,
        min_water_width: float = 0.1,
        roughness_mapping: List = ROUGHNESS_MAPPING_LIST,
    ) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame]:
        def rectangular_point_profile(
            branch: LineString,
            profiel_nummer: str,
            params: dict,
            def_height: float = 2,
            rect_offset: float = 0.1,
            interp_range: float = 0.1,
            thalweg_offset: float = None,
        ) -> List[Point]:

            bwidth = params["bodembreedte"]
            rotated_line = affinity.rotate(
                branch.geometry, 90, origin=branch.geometry.interpolate(0.5, normalized=True)
            )

            p1 = rotated_line.interpolate(0.5 - interp_range / 2, normalized=True)
            p2 = rotated_line.interpolate(0.5 + interp_range / 2, normalized=True)

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
                dx_1 = bwidth / np.sqrt(1 + s**2)
                dy_1 = dx_1 * s

                dx_2 = (bwidth + rect_offset) / np.sqrt(1 + s**2)
                dy_2 = dx_2 * s
            else:
                dx_1 = dx_2 = 0
                dy_1 = dy_2 = bwidth

            centroid = branch.geometry.interpolate(0.5, normalized=True)
            # centroid = rotated_line.centroid
            c_x, c_y = centroid.x, centroid.y

            if (thalweg_offset is not None) and (thalweg_offset > 0) and (thalweg_offset < 1):
                offset_l = thalweg_offset
                offset_r = 1 - thalweg_offset
            else:
                offset_l = offset_r = 0.5

            list_of_points = [
                Point(c_x - offset_l * dx_2, c_y - offset_l * dy_2),
                Point(c_x - offset_l * dx_1, c_y - offset_l * dy_1),
                Point(c_x + offset_r * dx_1, c_y + offset_r * dy_1),
                Point(c_x + offset_r * dx_2, c_y + offset_r * dy_2),
            ]

            bheight = np.nanmean(
                [params["bodemhoogte benedenstrooms"], params["bodemhoogte bovenstrooms"]]
            )
            if ("hoogte insteek linkerzijde" in params) and (
                "hoogte insteek rechterzijde" in params
            ):
                prof_height = [
                    params["hoogte insteek linkerzijde"],
                    bheight,
                    bheight,
                    params["hoogte insteek rechterzijde"],
                ]
            else:
                prof_height = [def_height, bheight, bheight, def_height]

            p_list = []
            for ix, point in enumerate(list_of_points):
                p_dict = dict(
                    [
                        ("codevolgnummer", ix + 1),
                        ("geometry", Point(point.x, point.y, prof_height[ix])),
                        ("profiel nummer", profiel_nummer),
                        ("type meting", 2),
                        ("typeruwheid", branch["typeruwheid"]),
                        ("ruwheidhoog", branch["ruwheidhoog"]),
                        ("ruwheidlaag", branch["ruwheidlaag"]),
                    ]
                )
                p_list.append(p_dict)

            return p_list

def load_geo_file(
    file_path: Any, geom_type: str = None, has_z_coord: bool = None, layer: str = None
) -> gpd.GeoDataFrame:
    """
    This function reads in a geospatial file and returns a geodataframe.

    When multiple paths are given, the files are combined.
    When given a list, a maximum of two files can be sjoined.
    When given a dict, multiple files can be sjoined or concatted, depending on the argument supplied in the dict.

    Args:
      file_path (Any): The path to the file you want to load. Can be multiple files (list or dict)
      layer (str): the name of the layer in the geopackage (if file_path ends in .gpkg). Defaults to None

    Errors:
      ValueError: when specific path does not end in either .shp or .gpkg, a ValueError is raised.
    Returns:
      gdf (gpd.GeoDataFrame): loaded ans possibly combined geodataframe
    """

    if isinstance(file_path, dict):
        for key, value in file_path.items():
            if "base" in key:
                gdf = load_geo_file(file_path=value)
            elif "concat" in key:
                gdf = gpd.GeoDataFrame(
                    data=pd.concat([gdf, load_geo_file(file_path=value)]),
                    geometry="geometry",
                    crs=gdf.crs,
                )
            elif "sjoin" in key:
                gdf2 = load_geo_file(file_path=value).to_crs(gdf.crs)
                gdf2.geometry = gdf2.geometry.buffer(5)
                gdf = gdf.sjoin(gdf2, how="left", predicate="within", rsuffix="right")
                gdf.columns = [column.replace("_right", "") for column in gdf.columns]

    elif isinstance(file_path, list):
        if len(file_path) <= 2:
            gdf_1 = load_geo_file(file_path=file_path[0])
            gdf_2 = load_geo_file(file_path=file_path[1])

            gdf_2_buffer = gdf_2.copy()
            gdf_2_buffer.geometry = gdf_2_buffer.geometry.buffer(1)
            gdf = gdf_1.sjoin(gdf_2_buffer, how="left", predicate="within")

        else:
            raise ValueError("Can only join two shapefiles together")

    elif file_path.endswith(r".shp"):
        gdf = gpd.read_file(file_path)
    elif file_path.endswith(r".gpkg"):
        if layer is not None:
            gdf = gpd.read_file(file_path, layer=layer)
        else:
            raise ValueError("provide a layer when loading a gpkg")
    else:
        raise ValueError("filetype not implemented")

    gdf.set_geometry(col="geometry", inplace=True)
    gdf = gdf.loc[~gdf["geometry"].isna(), :]

    # check for multi-element geometries
    if np.sum(gdf.geometry.type.isin(["MultiPoint", "MultiLine", "MultiLineString"])) > 0:
        gdf = gdf.explode(ignore_index=True, index_parts=False)

    if geom_type is not None:
        # space for additional checks
        if geom_type == "Point":
            if (gdf.geometry.type != "Point").any():
                pass
        elif geom_type == "LineString":
            if (gdf.geometry.type != "LineString").any():
                pass
        elif geom_type == "Polygon":
            if (gdf.geometry.type != "Polygon").any():
                pass
        else:
            raise ValueError("Unkown geometry type provided: {}".format(geom_type))

    if has_z_coord is not None:
        if not has_z_coord:
            if geom_type == "LineString":
                out_gdf = copy(gdf)
                for name, row in gdf.iterrows():
                    if np.array(row.geometry.coords).shape[1] > 2:
                        out_gdf.at[name, "geometry"] = LineString(
                            [xy[0:2] for xy in list(row.geometry.coords)]
                        )
                gdf = out_gdf
    return gdf.reset_index()


def map_columns(
    defaults, gdf: gpd.GeoDataFrame, index_mapping: dict, code_pad: str = None
) -> gpd.GeoDataFrame:
    """
    Function that:
        1. fixes possibly incorrect geometries in fix_geometry
        2. adds globalid columns. Overwrites possibly existing values to
        3. fills empty columns with values in fill_empty_columns
        4. checks for duplicate codes

    Args:
      defaults: a module of classes with default values for the columns in the GeoDataFrame
      gdf (gpd.GeoDataFrame): gpd.GeoDataFrame
      index_mapping (dict): dict
      code_pad (str): Padding used before a code. Defaults to None

    Returns:
      A GeoDataFrame with filled values and fixed geometries
    """

    def fill_empty_columns(
        defaults, gdf: gpd.GeoDataFrame, index_mapping: dict
    ) -> Tuple[gpd.GeoDataFrame, dict]:
        """
        Function that takes a GeoDataFrame, a dictionary with column names and a module of classes with default values. It then
        checks if the column names in the dictionary are present in the GeoDataFrame. If not, it adds
        the column with the default values. If the column is present, but contains missing values, it
        fills the missing values with the default values

        Args:
          defaults: a module with classes with default values for the columns that are not in the shapefile
          gdf (gpd.GeoDataFrame): gpd.GeoDataFrame
          index_mapping (dict): a dictionary with the DHYMAMO parameter name as key and the column name
        in the shapefile as value.

        Returns:
          A tuple of two objects: the tidy GeoDataFrame and a dictionary.
        """

        _index_mapping = copy(index_mapping)
        _gdf = copy(gdf)

        # loop over column names in index mapping, with DHYDAMO name in key and column name in shapefile in Value
        for key, value in index_mapping.items():
            _key = key.replace(" ", "_")
            # if data is not in shapefile, fill with default values
            if key == "globalid":
                continue
            elif value is None:
                if (defaults is not None) and hasattr(defaults, key):
                    _gdf[key] = getattr(defaults, _key)
                    _index_mapping[key] = key
                else:
                    _gdf[key] = None
                    _index_mapping[key] = None

            else:
                if isinstance(value, list):
                    try:
                        _gdf[key] = (
                            gdf[value]
                            .replace(to_replace=[None, "n.v.t.", -999, 0], value=np.nan)
                            .astype(float)
                            .apply(np.nanmean, axis=1)
                        )
                    except:
                        _gdf[key] = (
                            gdf[value]
                            .replace(to_replace=[None, "n.v.t.", -999, 0], value=np.nan)
                            .astype(str)
                            .bfill(axis=1)
                            .iloc[:, 0]  # Select the non-nan values
                        )
                elif isinstance(value, (int, float, bool)):
                    _gdf[key] = value
                else:
                    _gdf[key] = gdf[value]

                # if it is in shapefile, but contains missing values, fill them as well
                if key == "code":
                    duplicated = _gdf[key].duplicated(keep=False)
                    codes = [str(uuid.uuid4())[:8] for _ in range(np.sum(duplicated))]
                    _gdf.loc[duplicated, key] = codes

                    nones = _gdf[key].isna()
                    codes = [str(uuid.uuid4())[:8] for _ in range(np.sum(nones))]
                    _gdf.loc[nones, key] = codes

                elif _gdf[key].dtype == "geometry":
                    continue

                # elif key == "globalid":
                #     continue

                elif (_gdf[key].dtype == "int64") or (_gdf[key].dtype == "float64"):
                    if (defaults is not None) and hasattr(defaults, key):
                        _gdf[key] = _gdf[key].replace(
                            to_replace=[0, -999, np.nan],
                            value=getattr(defaults, _key),
                        )
                    # print number of 0, -998 and nan if no default value is present
                    else:
                        n_zero = _gdf[key][_gdf[key] == 0].count()
                        n_nnn = _gdf[key][_gdf[key] == -999].count()
                        n_nan = _gdf[key].isna().sum()

                        print(
                            "{}: entries with missing data: {}".format(key, n_zero + n_nnn + n_nan)
                        )

                elif (_gdf[key].dtype == "string") or (_gdf[key].dtype == "object"):
                    if hasattr(defaults, _key):
                        _gdf[key] = _gdf[key].replace(
                            to_replace=[None],
                            value=getattr(defaults, _key),
                        )
                    else:
                        n_nan = _gdf[key].isna().sum()
                        print("{}: entries with missing data: {}".format(key, n_nan))
                elif _gdf[key].dtype == bool:
                    pass
                else:
                    raise ValueError("type not yet implemented: {}".format(_gdf[key].dtype))

        return _gdf, _index_mapping

    def fix_geometry(gdf: gpd.GeoDataFrame, csr="epsg:28992") -> gpd.GeoDataFrame:
        """
        It drops rows without geometry and sets the correct crs

        Args:
          gdf (gpd.GeoDataFrame): GeoDataFrame to be fixed

        Returns:
          gdf (gpd.GeoDataFrame): GeoDataFrame with the geometry column dropped and the crs set to epsg:28992.
        """
        gdf.dropna(axis=0, inplace=True, subset=["geometry"])
        try:
            gdf = gdf.to_crs(csr)
        except ValueError:
            gdf = gdf.set_crs(csr).to_crs(csr)
        return gdf

    # ensure no incorrect geometries exist and set a globalid
    _gdf = fix_geometry(gdf)

    # fill potential empty columns
    _gdf, index_mapping = fill_empty_columns(
        defaults=defaults,
        gdf=_gdf,
        index_mapping=index_mapping,
    )
    _gdf["globalid"] = [str(uuid.uuid4()) for _ in range(_gdf.shape[0])]

    if ("code" in _gdf.columns) and (code_pad is not None):
        _gdf["code"] = code_pad + _gdf["code"].astype(str)
        # check for dupliacte columns in code
        _gdf = check_and_fix_duplicate_code(_gdf)

    # gdf = gpd.GeoDataFrame(_gdf[list(index_mapping.keys())], geometry="geometry", crs=gdf.crs)
    gdf = _gdf[list(index_mapping.keys())]

    # return only columns that are in that dict
    return gdf, index_mapping


def merge_multiple_to_dm(
    dm: Datamodel, features: List[str], gdfs: List[gpd.GeoDataFrame]
) -> Datamodel:
    if (not isinstance(features, List)) or (not isinstance(gdfs, List)):
        raise TypeError("inputs should be lists")

    for name, gdf in zip(features, gdfs):
        dm = merge_to_dm(dm, feature=name, feature_gdf=gdf)
    return dm


def merge_to_ddm(ddm: Datamodel, feature: str, feature_gdf: gpd.GeoDataFrame) -> Datamodel:
    """Legacy function, see merge_to_dm"""
    return merge_to_dm(dm=ddm, feature=feature, feature_gdf=feature_gdf)


def merge_to_dm(dm, feature: str, feature_gdf: gpd.GeoDataFrame):
    """
    Function that merges a feature into the supplied DataModel. If the feature already exists in the data model,
    then concatenate the new feature with the existing feature. Otherwise, add the new feature to the data model

    Args:
      dm (DataModel): the DataModel object object
      feature (str): name of the feature to be added to the DataModel
      feature_gdf (gpd.GeoDataFrame): the GeoDataFrame you want to merge into the DataModel

    Returns:
      dm (DataModel): A DataModel object with the new feature added.
    """
    if hasattr(dm, feature):
        new_gdf = gpd.GeoDataFrame(
            data=pd.concat([getattr(dm, feature), feature_gdf]),
            geometry="geometry",
            crs=feature_gdf.crs,
        )
        setattr(dm, feature, new_gdf)
    else:
        setattr(dm, feature, feature_gdf)

    return dm


def convert_to_dhydamo_data(dm: Datamodel, defaults: str, config: str) -> Datamodel:
    """
    This function creates a full DHYDAMO DataModel, based on an empty DataModel and a defaults and config file. It creates the
    DataModel by adding elements in the following order.
    - branches
    - culverts and bridges
    - measured profiles
    - peilgebieden or default peilen
    - pumps
    - river profiles
    - sluices and weirs

    Args:
      dm (Datamodel): empty DataModel, to be filled
      defaults (str): The path to the default file (should be in ./dataset_configs/)
      config (str): The path to the config file (should be in ./dataset_configs/).
    Returns:
      ddm (Datamodel): a fully HYDAMO compliant DataModel
    """

    def add_default_peil_to_branch(
        branches_gdf: gpd.GeoDataFrame, default_peil: float = None
    ) -> gpd.GeoDataFrame:
        """
        Add default peil values to branches.

        Args:
          branches_gdf (gpd.GeoDataFrame): a GeoDataFrame with a column called "peil"
          default_peil (float): value of default peil. Defaults to None

        Returns:
          A GeoDataFrame with a column "peil"
        """

        if default_peil is None:
            return branches_gdf

        if "peil" not in branches_gdf.columns.tolist():
            branches_gdf["peil"] = default_peil
        else:
            branches_gdf.loc[branches_gdf["peil"].isna(), "peil"] = default_peil
        return branches_gdf

    def add_peil_to_branch(
        branches_gdf: gpd.GeoDataFrame, peil_gdf: gpd.GeoDataFrame
    ) -> gpd.GeoDataFrame:
        in_cols = branches_gdf.columns.tolist()
        peil_gdf[["boven peil", "onder peil", "vast peil"]] = peil_gdf[["boven peil", "onder peil", "vast peil"]].replace(0, np.nan) #replace 0 by nan value
        peil_gdf["peil"] = peil_gdf[["boven peil", "onder peil", "vast peil"]].mean(axis=1, skipna=True) #peil is the average of bovenpeil, onderpeil and vastpeil

        #add ini peilen
        # peil_gdf["inipeil_zomer"]= None
        # peil_gdf["inipeil_winter"]= None      
        # for ix, waterloop in peil_gdf.iterrows():
        #     if ~ np.isnan(waterloop["vast peil"]):
        #         peil_gdf.loc[ix, "inipeil_zomer"] = waterloop["vast peil"]
        #         peil_gdf.loc[ix, "inipeil_winter"]= waterloop["vast peil"]
        #     elif ~ np.isnan(waterloop["onder peil"]) and ~ np.isnan(waterloop["boven peil"]):
        #         peil_gdf.loc[ix, "inipeil_zomer"] = waterloop["boven peil"]
        #         peil_gdf.loc[ix, "inipeil_winter"]= waterloop["onder peil"]
        #     else:
        #         print('Area has no peil: vast peil',waterloop["vast peil"],' onderpeil',waterloop["onder peil"]," bovenpeil",waterloop["boven peil"])
        #         #raise ValueError("Branches without peil")          

        old_geom = copy(branches_gdf["geometry"])
        branches_gdf["geometry"] = branches_gdf["geometry"].interpolate(0.5, normalized=True)
        out_branches = gpd.sjoin(
            left_df=branches_gdf,
            right_df=peil_gdf[["peil", "geometry"]].to_crs(branches_gdf.crs),
            how="left",
            predicate="within",
        )

        dups = out_branches.duplicated(subset="globalid", keep=False)
        if np.sum(dups) > 0:
            print(out_branches.loc[dups, :])
            raise ValueError("duplicates found")

        out_branches["geometry"] = old_geom
        in_cols.append("peil")
        #in_cols.append("inipeil_zomer")
        #in_cols.append("inipeil_winter")
        return out_branches[in_cols]

    def create_bridge_data(bridge_gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        """ """
        _bridge_gdf = copy(bridge_gdf)
        for ix, bridge in bridge_gdf.iterrows():
            if np.isnan(bridge["breedte_overspanning"]) or np.isnan(bridge["hoogte_onderkant"]):
                _bridge_gdf.drop(index=ix, inplace=True)
                continue
            # turn numerical roughnes types to strings
            _bridge_gdf.loc[ix, "typeruwheid"] = check_roughness(bridge)
            _bridge_gdf["doorstroomopening"] = None #van Laura

        for ix, bridge in _bridge_gdf.iterrows():
        #calculate the height of the cross section flow area
            hoogte_doorstroomopening = bridge.hoogte_onderkant - bridge["shift"] #both hoogte_onderkant and shift are wrt NAP
            crosssection = {
                    "shape": "rectangle",
                    "height": hoogte_doorstroomopening,
                    "width": bridge.breedte_overspanning,
                    "closed": 'no',
                }
            _bridge_gdf.loc[ix, "doorstroomopening"] = str(crosssection)

        return _bridge_gdf
    

    def create_culvert_data(culvert_gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        """ """

        culvert_gdf.set_index("code", inplace=True)
        culvert_gdf["doorstroomopening"] = None
        culvert_gdf["breedteopening"] = check_column_is_numerical(
            gdf=culvert_gdf["breedteopening"]
        )
        culvert_gdf["hoogteopening"] = check_column_is_numerical(gdf=culvert_gdf["hoogteopening"])
        culvert_gdf["lengte"] = check_column_is_numerical(gdf=culvert_gdf["lengte"])
        culvert_gdf["breedteopening"] = culvert_gdf["breedteopening"].replace(0, np.nan)
        culvert_gdf["hoogteopening"] = culvert_gdf["hoogteopening"].replace(0, np.nan)

        _culvert_gdf = copy(culvert_gdf)


        for ix, culvert in culvert_gdf.iterrows():
            if (not check_is_not_na_number(culvert["lengte"])) or (culvert["lengte"] == 1):
                if isinstance(culvert["geometry"], LineString):
                    _culvert_gdf.loc[ix, "lengte"] = culvert["geometry"].length

            if np.isnan(culvert["breedteopening"]) and np.isnan(culvert["hoogteopening"]):
                _culvert_gdf.drop(index=ix, inplace=True)
                continue
            elif np.isnan(culvert["breedteopening"]) or np.isnan(culvert["hoogteopening"]):
                shape = "circle"
                _culvert_gdf.loc[ix, "vormkoker"] = 1

            # elif (culvert["vormkoker"].dtype == "int64") or (culvert["vormkoker"].dtype == "float64"):
            elif isinstance(culvert["vormkoker"], int) | isinstance(culvert["vormkoker"], float):
                if int(culvert["vormkoker"]) == 3:
                    shape = "rectangle"
                    _culvert_gdf.loc[ix, "vormkoker"] = 3
                else:
                    shape = "circle"
                    _culvert_gdf.loc[ix, "vormkoker"] = 1

            elif isinstance(culvert["vormkoker"], str):  # dtype not a number, so assuming string
                if str(culvert["vormkoker"]).lower() == "rechthoekig":
                    shape = "rectangle"
                    _culvert_gdf.loc[ix, "vormkoker"] = 3
                else:
                    shape = "circle"
                    _culvert_gdf.loc[ix, "vormkoker"] = 1
            else:
                raise ValueError("wrong datatype: {}".format(type(culvert["vormkoker"])))

            if shape == "rectangle":
                crosssection = {
                    "shape": shape,
                    "height": culvert["hoogteopening"],
                    "width": culvert["breedteopening"],
                    "closed": culvert["gesloten"],
                }
            elif shape == "circle":
                diameter = np.nanmean(culvert[["breedteopening", "hoogteopening"]].values[:])
                _culvert_gdf.loc[ix, "breedteopening"] = diameter
                _culvert_gdf.loc[ix, "hoogteopening"] = diameter

                crosssection = {
                    "shape": shape,
                    "diameter": diameter,
                }
            _culvert_gdf.loc[ix, "doorstroomopening"] = str(crosssection)

            # turn numerical roughnes types to strings
            _culvert_gdf.loc[ix, "typeruwheid"] = check_roughness(culvert)

        _culvert_gdf = _culvert_gdf.drop(columns="gesloten")
        return _culvert_gdf.reset_index()

    def create_measured_profile_data(
        profile_points_gdf: gpd.GeoDataFrame,
        dist_tol: float = 0.25,
        roughness_mapping: List = ROUGHNESS_MAPPING_LIST,
    ) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame]:
        # Load shapefile with profile points & group them by metingprof attribute
        # profile_points = gpd.read_file(profile_points_path)

        profile_points_gdf["ruwheidhoog"] = check_column_is_numerical(
            profile_points_gdf["ruwheidhoog"]
        )
        profile_points_gdf["ruwheidlaag"] = check_column_is_numerical(
            profile_points_gdf["ruwheidlaag"]
        )

        grouped_points = profile_points_gdf.groupby(by="profiel nummer")
        # initialize empty lists
        profile_groups = []
        profile_lines = []
        profile_points = []
        roughnes_profiles = []

        # Loop over the grouped profile points to create a line from the points
        count = 0
        for name, group in grouped_points:
            # sort points based on codevolg number
            sorted_group = group.sort_values("codevolgnummer")

            # if both bottom and sludge measurements are present, select the first
            if (
                sorted_group["type meting"].isin([1]).any()
                and sorted_group["type meting"].isin([2]).any()
            ):
                sorted_group.drop(
                    sorted_group[sorted_group["type meting"] == 2].index, inplace=True
                )

            # skip 'line' if it only has one point
            if sorted_group.shape[0] < 4:
                continue

            # Initialize name and list of points
            profile_group_id = str(uuid.uuid4())
            profile_line_id = str(uuid.uuid4())
            #

            list_of_points = []

            # Add info to profile_group table
            _profile_group = dict(
                [
                    ("globalid", profile_group_id),
                    ("brugid", None),
                    ("stuwid", None),
                    ("geometry", None),
                ]
            )
            # _profile_group = dict([("globalid", profile_line_id), ("geometry", None)])
            profile_groups.append(_profile_group)

            code_volg_nr = 0
            # add points to line
            for ix, row in sorted_group.iterrows():
                if type(row.geometry) == Point:
                    l_points = [row.geometry]
                elif type(row.geometry) == MultiPoint:
                    l_points = row.geometry.geoms

                for point in l_points:
                    # Check if profile points are too close together, and skip if that is the case
                    # if ix < sorted_group.shape[0]:
                    if code_volg_nr > 0:
                        p_0 = list_of_points[code_volg_nr - 1]
                        p_1 = point
                        p_dist = p_0.distance(p_1)

                        if p_dist < dist_tol:
                            continue

                    # Because sorted_group has been sorted on codevolgnr, we can assume sequentiallity
                    code_volg_nr += 1

                    # append point to list for line generation
                    list_of_points.append(point)

                    # create entry for point shape and append
                    # point_id = "AGV_" + str(row["code"])
                    point_code = profile_line_id + r"_" + str(code_volg_nr)
                    point_id = str(uuid.uuid4())
                    _point = dict(
                        [
                            ("code", point_code),
                            ("globalid", point_id),
                            ("profiellijnid", profile_line_id),
                            # ("codevolgnummer", row["codevolgnu"]),
                            ("codevolgnummer", code_volg_nr),
                            ("geometry", point),
                        ]
                    )
                    profile_points.append(_point)

                    # create roughness table entry
                    _roughness_profile = dict(
                        [
                            ("code", point_code),
                            ("profielpuntid", point_id),
                            # ("typeruwheid", roughness_mapping[int(row["typeruwheid"]) - 1]),
                            ("typeruwheid", check_roughness(row)),
                            ("ruwheidhoog", float(row["ruwheidhoog"])),
                            ("ruwheidlaag", float(row["ruwheidlaag"])),
                            ("geometry", None),
                        ]
                    )
                    roughnes_profiles.append(_roughness_profile)

            # Convert points to line
            profile_line = LineString(list_of_points)

            # add line to list
            _profile_line = dict(
                [
                    ("globalid", profile_line_id),
                    ("profielgroepid", profile_group_id),
                    ("geometry", profile_line),
                ]
            )
            profile_lines.append(_profile_line)

        # Create GeoDataFrame from list of dicts
        profile_groups = gpd.GeoDataFrame(profile_groups, geometry="geometry", crs=28992)
        profile_lines = gpd.GeoDataFrame(profile_lines, geometry="geometry", crs=28992)
        profile_points = gpd.GeoDataFrame(profile_points, geometry="geometry", crs=28992)
        roughnes_profiles = gpd.GeoDataFrame(roughnes_profiles, geometry="geometry", crs=28992)

        return profile_groups, profile_lines, profile_points, roughnes_profiles

    def create_mp_from_np(
        branches_gdf: gpd.GeoDataFrame,
        dist_tol: float = 0.25,
        min_water_width: float = 0.1,
        roughness_mapping: List = ROUGHNESS_MAPPING_LIST,
    ) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame]:
        def rectangular_point_profile(
            branch: LineString,
            profiel_nummer: str,
            params: dict,
            def_height: float = 2,
            rect_offset: float = 0.1,
            interp_range: float = 0.1,
            thalweg_offset: float = None,
        ) -> List[Point]:

            bwidth = params["bodembreedte"]

            # l_offset_line_1 = branch.geometry.parallel_offset(
            #     distance=(bwidth + rect_offset) / 2, side="left"
            # )
            # l_offset_line_2 = branch.geometry.parallel_offset(distance=bwidth / 2, side="left")
            # r_offset_line_1 = branch.geometry.parallel_offset(
            #     distance=(bwidth + rect_offset) / 2, side="right"
            # )
            # r_offset_line_2 = branch.geometry.parallel_offset(distance=bwidth / 2, side="right")

            rotated_line = affinity.rotate(
                branch.geometry, 90, origin=branch.geometry.interpolate(0.5, normalized=True)
            )

            p1 = rotated_line.interpolate(0.5 - interp_range / 2, normalized=True)
            p2 = rotated_line.interpolate(0.5 + interp_range / 2, normalized=True)

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
                dx_1 = bwidth / np.sqrt(1 + s**2)
                dy_1 = dx_1 * s

                dx_2 = (bwidth + rect_offset) / np.sqrt(1 + s**2)
                dy_2 = dx_2 * s
            else:
                dx_1 = dx_2 = 0
                dy_1 = dy_2 = bwidth

            centroid = branch.geometry.interpolate(0.5, normalized=True)
            # centroid = rotated_line.centroid
            c_x, c_y = centroid.x, centroid.y

            if (thalweg_offset is not None) and (thalweg_offset > 0) and (thalweg_offset < 1):
                offset_l = thalweg_offset
                offset_r = 1 - thalweg_offset
            else:
                offset_l = offset_r = 0.5

            list_of_points = [
                Point(c_x - offset_l * dx_2, c_y - offset_l * dy_2),
                Point(c_x - offset_l * dx_1, c_y - offset_l * dy_1),
                Point(c_x + offset_r * dx_1, c_y + offset_r * dy_1),
                Point(c_x + offset_r * dx_2, c_y + offset_r * dy_2),
            ]

            bheight = np.nanmean(
                [params["bodemhoogte benedenstrooms"], params["bodemhoogte bovenstrooms"]]
            )
            if ("hoogte insteek linkerzijde" in params) and (
                "hoogte insteek rechterzijde" in params
            ):
                prof_height = [
                    params["hoogte insteek linkerzijde"],
                    bheight,
                    bheight,
                    params["hoogte insteek rechterzijde"],
                ]
            else:
                prof_height = [def_height, bheight, bheight, def_height]

            p_list = []
            for ix, point in enumerate(list_of_points):
                p_dict = dict(
                    [
                        ("codevolgnummer", ix + 1),
                        ("geometry", Point(point.x, point.y, prof_height[ix])),
                        ("profiel nummer", profiel_nummer),
                        ("type meting", 2),
                        ("typeruwheid", branch["typeruwheid"]),
                        ("ruwheidhoog", branch["ruwheidhoog"]),
                        ("ruwheidlaag", branch["ruwheidlaag"]),
                    ]
                )
                p_list.append(p_dict)

            return p_list

        def trapezium_point_profile(
            branch: LineString,
            profiel_nummer: str,
            params: dict,
            interp_range: float = 0.1,
            thalweg_offset: float = None,
        ) -> List[Point]:

            bheight = np.nanmean(
                [params["bodemhoogte benedenstrooms"], params["bodemhoogte bovenstrooms"]]
            )

            prof_height = [
                params["hoogte insteek linkerzijde"],
                bheight,
                bheight,
                params["hoogte insteek rechterzijde"],
            ]

            l_depth = params["hoogte insteek linkerzijde"] - bheight
            r_depth = params["hoogte insteek rechterzijde"] - bheight

            # talud is dx/dy
            left_offset = params["taludhelling linkerzijde"] * l_depth
            right_offset = params["taludhelling rechterzijde"] * r_depth

            bwidth = params["bodembreedte"]

            # l_offset_line_1 = branch.geometry.parallel_offset(
            #     distance=bwidth / 2 + left_offset, side="left"
            # )
            # l_offset_line_2 = branch.geometry.parallel_offset(distance=bwidth / 2, side="left")
            # r_offset_line_1 = branch.geometry.parallel_offset(
            #     distance=bwidth / 2 + right_offset, side="right"
            # )
            # r_offset_line_2 = branch.geometry.parallel_offset(distance=bwidth / 2, side="right")

            rotated_line = affinity.rotate(
                branch.geometry, 90, origin=branch.geometry.interpolate(0.5, normalized=True)
            )

            p1 = rotated_line.interpolate(0.5 - interp_range / 2, normalized=True)
            p2 = rotated_line.interpolate(0.5 + interp_range / 2, normalized=True)

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
                dx_1 = bwidth / np.sqrt(1 + s**2)
                dy_1 = dx_1 * s

                dx_2 = (bwidth + left_offset + right_offset) / np.sqrt(1 + s**2)
                dy_2 = dx_2 * s
            else:
                dx_1 = dx_2 = 0
                dy_1 = dy_2 = bwidth

            centroid = branch.geometry.interpolate(0.5, normalized=True)
            # centroid = rotated_line.centroid
            c_x, c_y = centroid.x, centroid.y

            if (thalweg_offset is not None) and (thalweg_offset > 0) and (thalweg_offset < 1):
                offset_l = thalweg_offset
                offset_r = 1 - thalweg_offset
            else:
                offset_l = offset_r = 0.5

            list_of_points = [
                Point(c_x - offset_l * dx_2, c_y - offset_l * dy_2),
                Point(c_x - offset_l * dx_1, c_y - offset_l * dy_1),
                Point(c_x + offset_r * dx_1, c_y + offset_r * dy_1),
                Point(c_x + offset_r * dx_2, c_y + offset_r * dy_2),
            ]

            p_list = []
            for ix, point in enumerate(list_of_points):
                try:
                    p_dict = dict(
                        [
                            ("codevolgnummer", ix + 1),
                            ("geometry", Point(point.x, point.y, prof_height[ix])),
                            ("profiel nummer", profiel_nummer),
                            ("type meting", 2),
                            ("typeruwheid", branch["typeruwheid"]),
                            ("ruwheidhoog", branch["ruwheidhoog"]),
                            ("ruwheidlaag", branch["ruwheidlaag"]),
                        ]
                    )
                except:
                    print(params)
                    raise
                p_list.append(p_dict)
            return p_list

        # Load shapefile with profile points & group them by metingprof attribute
        # profile_points = gpd.read_file(profile_points_path)
        branches_gdf[
            [
                "bodembreedte",
                "bodemhoogte benedenstrooms",
                "bodemhoogte bovenstrooms",
                "hoogte insteek linkerzijde",
                "hoogte insteek rechterzijde",
                "taludhelling linkerzijde",
                "taludhelling rechterzijde",
            ]
        ] = branches_gdf[
            [
                "bodembreedte",
                "bodemhoogte benedenstrooms",
                "bodemhoogte bovenstrooms",
                "hoogte insteek linkerzijde",
                "hoogte insteek rechterzijde",
                "taludhelling linkerzijde",
                "taludhelling rechterzijde",
            ]
        ].astype(
            float
        )

        branches_gdf.loc[
            branches_gdf["bodembreedte"] < min_water_width,
            "bodembreedte",
        ] = np.nan
        branches_gdf.loc[
            branches_gdf["water_width_index"] < min_water_width,
            "water_width_index",
        ] = np.nan

        branches_gdf["globalid"] = [str(uuid.uuid4()) for _ in range(branches_gdf.shape[0])]
        branches_gdf.set_index("globalid", drop=False, inplace=True)
        branches_out_gdf = copy(branches_gdf)

        # initialize empty lists
        pp_list = []
        # Loop over the grouped profile points to create a line from the points
        for jx, ix_1 in enumerate(branches_gdf.index):
            branch = copy(branches_gdf.loc[ix_1, :])
            # Dont use jx with numerical index column

            # # Create unique ident for branch and add to branch
            branch_gid = ix_1

            # # replace duplicate codes by addinng a uniqure number
            # if duplicates.iloc[jx]:
            #     branches_out_gdf.at[ix_1, "code"] = branches_out_gdf.at[ix_1, "code"] + str(jx)

            # turn numerical roughnes types to strings
            type_ruwheid = check_roughness(branch)
            branches_out_gdf.at[ix_1, "typeruwheid"] = type_ruwheid

            # check if linestring
            try:
                assert isinstance(branch.geometry, LineString)
            except AssertionError:
                print(branch.geometry.type)
                raise

            # Delete z-dimensison of linestring if it exists
            if np.array(branch.geometry.coords).shape[1] > 2:
                branches_out_gdf.at[ix_1, "geometry"] = LineString(
                    [xy[0:2] for xy in list(branch.geometry.coords)]
                )

            # Check if width and depth parameters are available, if not, skip
            # Also skip if no width is available
            if (
                (
                    (branch[["bodemhoogte benedenstrooms", "bodemhoogte bovenstrooms"]].empty)
                    or (
                        branch[["bodemhoogte benedenstrooms", "bodemhoogte bovenstrooms"]]
                        .isna()
                        .values.all()
                    )
                )
                # or (
                #     (
                #         (branch[["bodembreedte", "bodemhoogte benedenstrooms"]].empty)
                #         or (
                #             branch[["bodembreedte", "bodemhoogte benedenstrooms"]]
                #             .isna()
                #             .values.all()
                #         )
                #     )
                #     and (
                #         (branch[["water_width_index", "bodemhoogte benedenstrooms"]].empty)
                #         or (
                #             branch[["water_width_index", "bodemhoogte benedenstrooms"]]
                #             .isna()
                #             .values.all()
                #         )
                #     )
                # )
                or (
                    (branch[["bodembreedte", "water_width_index"]].empty)
                    or (branch[["bodembreedte", "water_width_index"]].isna().values.all())
                )
            ):
                # print("no profile for branch {}".format(branch_gid))
                # print(branch)
                continue

            if (
                (not branch[["bodembreedte", "water_width_index"]].isna().values.any())
                and (not branch[["bodembreedte"]].empty)
                and (not branch[["water_width_index"]].empty)
            ):
                if branch["bodembreedte"] > branch["water_width_index"]:
                    branch["bodembreedte"] = np.nan

            # Check if all parameters are present for trapezium profile. If not, use rectangular profile
            width_ix = "bodembreedte"
            if np.isnan(branch["bodembreedte"]):
                prof_type = "rectangle"
                width_ix = "water_width_index"
            elif (
                (not check_is_not_na_number(branch["hoogte insteek linkerzijde"]))
                | (not check_is_not_na_number(branch["hoogte insteek rechterzijde"]))
                | (not check_is_not_na_number(branch["taludhelling linkerzijde"]))
                | (not check_is_not_na_number(branch["taludhelling rechterzijde"]))
            ):
                prof_type = "rectangle"
            else:
                prof_type = "trapezium"

            width = branch[width_ix]
            if width < 0:
                width = 0
            # elif width > 200:
            #     width = 0

            # set required parameters for either rectangle or trapezium profile
            if prof_type == "rectangle":
                params = dict(
                    [
                        ("bodembreedte", width),
                        (
                            "bodemhoogte benedenstrooms",
                            branch["bodemhoogte benedenstrooms"],
                        ),
                        (
                            "bodemhoogte bovenstrooms",
                            branch["bodemhoogte bovenstrooms"],
                        ),
                    ]
                )
                if check_is_not_na_number(branch["hoogte insteek linkerzijde"]):
                    params["hoogte insteek linkerzijde"] = branch["hoogte insteek linkerzijde"]

                if check_is_not_na_number(branch["hoogte insteek rechterzijde"]):
                    params["hoogte insteek rechterzijde"] = branch["hoogte insteek rechterzijde"]

                if ("thalweg offset" in branch) and check_is_not_na_number(
                    branch["thalweg offset"]
                ):
                    thalweg_offset = branch["thalweg offset"]
                else:
                    thalweg_offset = None

                p_list = rectangular_point_profile(
                    branch=branch,
                    profiel_nummer=ix_1,
                    params=params,
                    thalweg_offset=thalweg_offset,
                )
                pp_list += p_list
            elif prof_type == "trapezium":
                params = dict(
                    [
                        ("bodembreedte", width),
                        (
                            "bodemhoogte benedenstrooms",
                            branch["bodemhoogte benedenstrooms"],
                        ),
                        (
                            "bodemhoogte bovenstrooms",
                            branch["bodemhoogte bovenstrooms"],
                        ),
                        (
                            "hoogte insteek linkerzijde",
                            branch["hoogte insteek linkerzijde"],
                        ),
                        (
                            "hoogte insteek rechterzijde",
                            branch["hoogte insteek rechterzijde"],
                        ),
                        (
                            "taludhelling linkerzijde",
                            branch["taludhelling linkerzijde"],
                        ),
                        (
                            "taludhelling rechterzijde",
                            branch["taludhelling rechterzijde"],
                        ),
                    ]
                )

                # if hoogte insteek is lower than the bottom, swap bottom and hoogte insteek
                if np.amin(
                    [params["hoogte insteek linkerzijde"], params["hoogte insteek rechterzijde"]]
                ) < np.amax(
                    [params["bodemhoogte benedenstrooms"], params["bodemhoogte bovenstrooms"]]
                ):
                    (
                        params["bodemhoogte benedenstrooms"],
                        params["bodemhoogte bovenstrooms"],
                        params["hoogte insteek linkerzijde"],
                        params["hoogte insteek rechterzijde"],
                    ) = (
                        params["hoogte insteek linkerzijde"],
                        params["hoogte insteek rechterzijde"],
                        params["bodemhoogte benedenstrooms"],
                        params["bodemhoogte bovenstrooms"],
                    )

                if ("thalweg offset" in branch) and check_is_not_na_number(
                    branch["thalweg offset"]
                ):
                    thalweg_offset = branch["thalweg offset"]
                else:
                    thalweg_offset = None

                p_list = trapezium_point_profile(
                    branch=branch,
                    profiel_nummer=ix_1,
                    params=params,
                    thalweg_offset=thalweg_offset,
                )
                pp_list += p_list

        # Create GeoDataFrame from list of dicts
        pp_gdf = gpd.GeoDataFrame(data=pp_list, geometry="geometry", crs=branches_gdf.crs)
        return create_measured_profile_data(profile_points_gdf=pp_gdf, dist_tol=0)

    def create_norm_parm_profiles_v2(
        branches_gdf: gpd.GeoDataFrame, min_water_width: float = 0.1
    ) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame]:
        """
        Function that converts parameterized profiles to HYDAMO compliant format.

        Args:
            branches_gdf (gpd.GeoDataFrame): input geodataframe containing branches
            index_mapping (dict): dictionary containing a mapping from required keys to values present in branches_gdf

        Returns:
            hydroobject (gpd.GeoDataFrame): ouput geodataframe containing branches
            hydroobject_normgp (gpd.GeoDataFrame): output geodataframe containing names of branches and corresponding profiles
            normgeparamprofielwaarde (gpd.GeoDataFrame): output geodataframe containing parameterized profiles
        """
        # script
        ho_ngp_list = []
        ngp_list = []
        branches_gdf[
            [
                "bodembreedte",
                "bodemhoogte benedenstrooms",
                "bodemhoogte bovenstrooms",
                "hoogte insteek linkerzijde",
                "hoogte insteek rechterzijde",
                "taludhelling linkerzijde",
                "taludhelling rechterzijde",
            ]
        ] = branches_gdf[
            [
                "bodembreedte",
                "bodemhoogte benedenstrooms",
                "bodemhoogte bovenstrooms",
                "hoogte insteek linkerzijde",
                "hoogte insteek rechterzijde",
                "taludhelling linkerzijde",
                "taludhelling rechterzijde",
            ]
        ].astype(
            float
        )

        branches_gdf.loc[
            branches_gdf["bodembreedte"] < min_water_width,
            "bodembreedte",
        ] = np.nan
        branches_gdf.loc[
            branches_gdf["water_width_index"] < min_water_width,
            "water_width_index",
        ] = np.nan

        branches_gdf["globalid"] = [str(uuid.uuid4()) for _ in range(branches_gdf.shape[0])]
        branches_gdf.set_index("globalid", drop=False, inplace=True)
        branches_out_gdf = copy(branches_gdf)

        # # Check for duplicate codes, if there are, they are replaced in the following loop
        # duplicates = branches_gdf.duplicated(subset=["code"])
        # for ix_1, branch in tqdm(branches_gdf.iterrows(), total=branches_gdf.shape[0]):
        # for jx, (ix_1, branch) in enumerate(branches_gdf.iterrows()):
        for jx, ix_1 in enumerate(branches_gdf.index):
            branch = copy(branches_gdf.loc[ix_1, :])
            # Dont use jx with numerical index column

            # # Create unique ident for branch and add to branch
            branch_gid = ix_1

            # # replace duplicate codes by addinng a uniqure number
            # if duplicates.iloc[jx]:
            #     branches_out_gdf.at[ix_1, "code"] = branches_out_gdf.at[ix_1, "code"] + str(jx)

            # turn numerical roughnes types to strings
            type_ruwheid = check_roughness(branch)
            branches_out_gdf.at[ix_1, "typeruwheid"] = type_ruwheid

            # check if linestring
            assert isinstance(branch.geometry, LineString)

            # Delete z-dimensison of linestring if it exists
            if np.array(branch.geometry.coords).shape[1] > 2:
                branches_out_gdf.at[ix_1, "geometry"] = LineString(
                    [xy[0:2] for xy in list(branch.geometry.coords)]
                )

            # Check if width and depth parameters are available, if not, skip
            # Also skip if no width is available
            if (
                (
                    (branch[["bodemhoogte benedenstrooms", "bodemhoogte bovenstrooms"]].empty)
                    or (
                        branch[["bodemhoogte benedenstrooms", "bodemhoogte bovenstrooms"]]
                        .isna()
                        .values.all()
                    )
                )
                # or (
                #     (
                #         (branch[["bodembreedte", "bodemhoogte benedenstrooms"]].empty)
                #         or (
                #             branch[["bodembreedte", "bodemhoogte benedenstrooms"]]
                #             .isna()
                #             .values.all()
                #         )
                #     )
                #     and (
                #         (branch[["water_width_index", "bodemhoogte benedenstrooms"]].empty)
                #         or (
                #             branch[["water_width_index", "bodemhoogte benedenstrooms"]]
                #             .isna()
                #             .values.all()
                #         )
                #     )
                # )
                or (
                    (branch[["bodembreedte", "water_width_index"]].empty)
                    or (branch[["bodembreedte", "water_width_index"]].isna().values.all())
                )
            ):
                # print("no profile for branch {}".format(branch_gid))
                continue

            if (
                (not branch[["bodembreedte", "water_width_index"]].isna().values.any())
                and (not branch[["bodembreedte"]].empty)
                and (not branch[["water_width_index"]].empty)
            ):
                if branch["bodembreedte"] > branch["water_width_index"]:
                    branch["bodembreedte"] = np.nan

            # Check if all parameters are present for trapezium profile. If not, use rectangular profile
            width_ix = "bodembreedte"
            if np.isnan(branch["bodembreedte"]):
                prof_type = "rectangle"
                width_ix = "water_width_index"
            # elif (
            #     (branch["hoogte insteek linkerzijde"] is None)
            #     | (branch["hoogte insteek rechterzijde"] is None)
            #     | (branch["taludhelling linkerzijde"] is None)
            #     | (branch["taludhelling rechterzijde"] is None)
            # ):
            #     prof_type = "rectangle"
            # elif (
            #     np.isnan(branch["hoogte insteek linkerzijde"])
            #     | np.isnan(branch["hoogte insteek rechterzijde"])
            #     | np.isnan(branch["taludhelling linkerzijde"])
            #     | np.isnan(branch["taludhelling rechterzijde"])
            # ):
            #     prof_type = "rectangle"
            # elif (
            #     (branch["hoogte insteek linkerzijde"] == 0)
            #     | (branch["hoogte insteek rechterzijde"] == 0)
            #     | (branch["taludhelling linkerzijde"] == 0)
            #     | (branch["taludhelling rechterzijde"] == 0)
            # ):
            elif (
                (not check_is_not_na_number(branch["hoogte insteek linkerzijde"]))
                | (not check_is_not_na_number(branch["hoogte insteek rechterzijde"]))
                | (not check_is_not_na_number(branch["taludhelling linkerzijde"]))
                | (not check_is_not_na_number(branch["taludhelling rechterzijde"]))
            ):
                prof_type = "rectangle"
            else:
                prof_type = "trapezium"

            width = branch[width_ix]
            if width < 0:
                width = 0
            # elif width > 200:
            #     width = 0

            # add entry to hydroobject_normgp tabel
            ngp_gid = str(uuid.uuid4())
            ngp = dict(
                [
                    ("hydroobjectid", branch_gid),
                    ("normgeparamprofielid", ngp_gid),
                    ("globalid", ngp_gid),
                    ("geometry", None),
                ]
            )
            ho_ngp_list.append(ngp)

            # set required parameters for either rectangle or trapezium profile
            if prof_type == "rectangle":
                params = dict(
                    [
                        ("bodembreedte", width),
                        (
                            "bodemhoogte benedenstrooms",
                            branch["bodemhoogte benedenstrooms"],
                        ),
                        (
                            "bodemhoogte bovenstrooms",
                            branch["bodemhoogte bovenstrooms"],
                        ),
                    ]
                )
            elif prof_type == "trapezium":
                params = dict(
                    [
                        ("bodembreedte", width),
                        (
                            "bodemhoogte benedenstrooms",
                            branch["bodemhoogte benedenstrooms"],
                        ),
                        (
                            "bodemhoogte bovenstrooms",
                            branch["bodemhoogte bovenstrooms"],
                        ),
                        (
                            "hoogte insteek linkerzijde",
                            branch["hoogte insteek linkerzijde"],
                        ),
                        (
                            "hoogte insteek rechterzijde",
                            branch["hoogte insteek rechterzijde"],
                        ),
                        (
                            "taludhelling linkerzijde",
                            branch["taludhelling linkerzijde"],
                        ),
                        (
                            "taludhelling rechterzijde",
                            branch["taludhelling rechterzijde"],
                        ),
                    ]
                )

                # if hoogte insteek is lower than the bottom, swap bottom and hoogte insteek
                if np.amin(
                    [params["hoogte insteek linkerzijde"], params["hoogte insteek rechterzijde"]]
                ) < np.amax(
                    [params["bodemhoogte benedenstrooms"], params["bodemhoogte bovenstrooms"]]
                ):
                    (
                        params["bodemhoogte benedenstrooms"],
                        params["bodemhoogte bovenstrooms"],
                        params["hoogte insteek linkerzijde"],
                        params["hoogte insteek rechterzijde"],
                    ) = (
                        params["hoogte insteek linkerzijde"],
                        params["hoogte insteek rechterzijde"],
                        params["bodemhoogte benedenstrooms"],
                        params["bodemhoogte bovenstrooms"],
                    )

            # loop over parameters to add to ngp_list
            for ix_2, (key, value) in enumerate(params.items()):
                ngp_values = dict(
                    [
                        ("normgeparamprofielid", ngp_gid),
                        (
                            "typeruwheid",
                            type_ruwheid,
                        ),
                        ("ruwheidhoog", branch["ruwheidhoog"]),
                        ("ruwheidlaag", branch["ruwheidlaag"]),
                        ("soortparameter", key),
                        ("waarde", value),
                        ("geometry", None),
                    ]
                )
                ngp_list.append(ngp_values)

        if len(ho_ngp_list) > 0:
            hydroobject_normgp = gpd.GeoDataFrame(
                ho_ngp_list, geometry="geometry", crs=branches_out_gdf.crs
            )
        else:
            hydroobject_normgp = None

        if len(ngp_list) > 0:
            normgeparamprofielwaarde = gpd.GeoDataFrame(
                ngp_list, geometry="geometry", crs=branches_out_gdf.crs
            )
        else:
            normgeparamprofielwaarde = None

        return (
            # branches_out_gdf[["code", "globalid", "geometry", "typeruwheid"]].reset_index(
            #     drop=True
            # ),
            branches_out_gdf.reset_index(drop=True),
            hydroobject_normgp,
            normgeparamprofielwaarde,
        )

    def create_pump_data(pump_gdf: gpd.GeoDataFrame,peil_gdf: gpd.GeoDataFrame) -> List[gpd.GeoDataFrame]:
        pump_station_list = []
        pump_list = []
        management_list = []

        pump_gdf["maximalecapaciteit"] = check_column_is_numerical(
            gdf=pump_gdf["maximalecapaciteit"]
        )

        pump_gdf=create_streefpeil(peil_gdf,pump_gdf)


        for ix, pump in pump_gdf.iterrows():
            pump_station = dict(
                [
                    ("code", pump["code"]),
                    ("geometry", pump["geometry"]),
                    ("globalid", pump["globalid"]),
                ]
            )
            pump_station_list.append(pump_station)

            pump_gid = str(uuid.uuid4())
            _pump = dict(
                [
                    ("code", pump["code"]),
                    ("geometry", None),
                    ("gemaalid", pump["globalid"]),
                    ("globalid", pump_gid),
                    ("maximalecapaciteit", pump["maximalecapaciteit"]),
                ]
            )
            pump_list.append(_pump)
     
            management = dict(
                [
                    ("bovengrens", pump["streefwaarde"] + 0.5 * pump["peil_marge"]),
                    ("code", pump["code"]),
                    ("doelvariabele", pump["doelvariabele"]),
                    ("geometry", None),
                    ("globalid", str(uuid.uuid4())),
                    ("ondergrens", pump["streefwaarde"] - 0.5 * pump["peil_marge"]),
                    ("pompid", pump_gid),
                    ("streefwaarde", pump["streefwaarde"]),
                ]
            )
            management_list.append(management)

        pump_station_gdf = gpd.GeoDataFrame(pump_station_list, geometry="geometry", crs=28992)
        _pump_gdf = gpd.GeoDataFrame(pump_list, geometry="geometry", crs=28992)
        management_gdf = gpd.GeoDataFrame(management_list, geometry="geometry", crs=28992)

        return pump_station_gdf, _pump_gdf, management_gdf

    def create_river_profiles(
        branches_gdf: gpd.GeoDataFrame,
        riv_prof_df: pd.DataFrame,
        code_padding: str,
        rough_df: pd.DataFrame = None,
    ) -> gpd.GeoDataFrame:

        riv_prof_df["id"] = riv_prof_df["id"].astype(str)
        unique_ids = riv_prof_df["id"].unique()

        geom_list = []
        meta_list = []

        for u_id in unique_ids:
            slice = riv_prof_df.loc[riv_prof_df["id"] == u_id, :]
            meta_data = slice.loc[slice["Data_type"] == "meta", :]
            geom_data = slice.loc[slice["Data_type"] == "geom", :]

            branch = str(meta_data["branch"].values[0])
            _branchid = code_padding + branch
            name = code_padding + u_id

            if np.sum(branches_gdf["code"].str.contains(_branchid)) == 0:
                continue

            total_width = (
                meta_data["width main channel"].values[0]
                + meta_data["width floodplain 1"].values[0]
                + meta_data["width floodplain 2"].values[0]
            )
            max_width = geom_data["Flow width"].values[0].max()
            if total_width > max_width:
                diff = total_width - max_width
            else:
                diff = None

            for ix in range(geom_data.shape[0]):
                flow_width = geom_data["Flow width"].values[ix]
                total_width = geom_data["Total width"].values[ix]

                if (flow_width == max_width) and (diff is not None):
                    flow_width = max_width + diff
                    total_width = max_width + diff + 0.1

                geom_dict = dict(
                    [
                        ("name", name),
                        ("ix", ix),
                        ("levels", geom_data["level"].values[ix]),
                        ("flowWidths", flow_width),
                        ("totalWidths", total_width),
                        ("geometry", None),
                    ]
                )
                geom_list.append(geom_dict)

            branchid = branches_gdf.loc[
                branches_gdf["code"].str.contains(_branchid),
                "code",
            ].values[0]

            meta_dict = dict(
                [
                    ("numLevels", geom_data.shape[0]),
                    ("name", name),
                    ("thalweg", 0),
                    ("leveecrestLevel", meta_data["Crest level summerdike"].values[0]),
                    (
                        "leveebaselevel",
                        meta_data["Floodplain baselevel behind summerdike"].values[0],
                    ),
                    ("leveeflowarea", meta_data["Flow area behind summerdike"].values[0]),
                    ("leveetotalarea", meta_data["Total area behind summerdike"].values[0]),
                    ("mainwidth", meta_data["width main channel"].values[0]),
                    ("fp1width", meta_data["width floodplain 1"].values[0]),
                    ("fp2width", meta_data["width floodplain 2"].values[0]),
                    ("branchid", branchid),
                    ("chainage", meta_data["chainage"].values[0]),
                    ("geometry", None),
                ]
            )
            if rough_df is not None:
                meta_dict["frictionids"] = "Main,FloodPlain1,FloodPlain1"
            else:
                meta_dict["frictionids"] = None

            meta_list.append(meta_dict)

        geom_gdf = gpd.GeoDataFrame(geom_list, geometry="geometry", crs="epsg:28992")
        meta_gdf = gpd.GeoDataFrame(meta_list, geometry="geometry", crs="epsg:28992")

        if rough_df is not None:
            rough_list = []
            branches = riv_prof_df.loc[riv_prof_df["Data_type"] == "meta", "branch"].unique()
            for ix, branch in enumerate(branches):
                if np.sum(branches_gdf["code"].str.contains(branch)) == 0:
                    continue

                rough_slice = rough_df.loc[rough_df["Name"] == branch, :]
                sections = ["Main", "FloodPlain1", "FloodPlain2"]
                branchid = branches_gdf.loc[
                    branches_gdf["code"].str.contains(branch),
                    "code",
                ].values[0]

                for jx, section in enumerate(sections):
                    section_slice = rough_slice.loc[rough_slice["SectionType"] == section, :]

                    if section_slice.empty:
                        continue

                    function_type = BRANCH_FRICTION_FUNCTION[section_slice["Dependance"].values[0]]
                    values = ",".join(map(str, section_slice["R_pos_constant"].values[:]))
                    slice_dict = dict(
                        [
                            ("branchid", branchid),
                            ("frictiontype", section_slice["RoughnessType"].values[0]),
                            ("functiontype", function_type),
                            ("geometry", None),
                            ("section", section),
                            ("frictionvalues", values),
                        ]
                    )

                    chainages = section_slice["Chainage"]
                    n_levels = chainages.value_counts().max()
                    locations = chainages.value_counts()

                    if n_levels > 1:
                        chainage_slice = section_slice.loc[
                            section_slice["Chainage"] == chainages.iloc[0], :
                        ]

                        slice_dict["numlevels"] = n_levels

                        if function_type == "absDischarge":
                            levels = chainage_slice["Q_pos"].values[:]
                            values = section_slice["R_pos_f(Q)"].values[:]
                            slice_dict["levels"] = ",".join(map(str, levels))
                            slice_dict["frictionvalues"] = ",".join(map(str, values))
                        elif function_type == "waterLevel":
                            levels = chainage_slice["H_pos"].values[:]
                            values = section_slice["R_pos__f(h)"].values[:]
                            slice_dict["levels"] = ",".join(map(str, levels))
                            slice_dict["frictionvalues"] = ",".join(map(str, values))
                        else:
                            print(function_type)
                            raise ValueError("unknown function type")

                    if locations.shape[0] > 1:
                        slice_dict["numlocations"] = locations.shape[0]
                        slice_dict["chainage"] = ",".join(map(str, locations.index))

                    rough_list.append(slice_dict)

            rough_gdf = gpd.GeoDataFrame(rough_list, geometry="geometry", crs="epsg:28992")

        else:
            rough_gdf = None

        return geom_gdf, meta_gdf, rough_gdf

    def create_weir_data(
        branches_gdf: gpd.GeoDataFrame, weir_gdf: gpd.GeoDataFrame
    ) -> List[gpd.GeoDataFrame]:
        if branches_gdf is None:
            branch_geom = None
        else:
            branch_geom = branches_gdf["geometry"]

        weir_list = []
        opening_list = []
        management_device_list = []

        # if not "hoogstedoorstroombreedte" in weir_gdf.columns:
        #     weir_gdf["hoogstedoorstroombreedte"] = weir_gdf["laagstedoorstroombreedte"]

        # if not "hoogstedoorstroomhoogte" in weir_gdf.columns:
        #     weir_gdf["hoogstedoorstroomhoogte"] = weir_gdf["laagstedoorstroomhoogte"] + 10
        clist = [
            "hoogstedoorstroombreedte",
            "hoogstedoorstroomhoogte",
            "laagstedoorstroombreedte",
            "laagstedoorstroomhoogte",
        ]
        for column in clist:
            weir_gdf[column] = check_column_is_numerical(gdf=weir_gdf[column])

        weir_gdf["hoogstedoorstroombreedte"].fillna(
            weir_gdf["laagstedoorstroombreedte"], inplace=True
        )
        weir_gdf["hoogstedoorstroomhoogte"].fillna(
            weir_gdf["laagstedoorstroomhoogte"], inplace=True
        )
        for ix, weir in weir_gdf.iterrows():
            if np.isnan(weir["laagstedoorstroombreedte"]) and np.isnan(
                weir["hoogstedoorstroombreedte"]
            ):
                continue
            elif np.isnan(weir["laagstedoorstroombreedte"]) or np.isnan(
                weir["hoogstedoorstroombreedte"]
            ):
                b_laag = b_hoog = np.nanmean(
                    weir[["hoogstedoorstroombreedte", "laagstedoorstroombreedte"]].values[:]
                )
            else:
                b_laag = weir["laagstedoorstroombreedte"]
                b_hoog = weir["hoogstedoorstroombreedte"]
            if np.isnan(weir["laagstedoorstroomhoogte"]) and np.isnan(
                weir["hoogstedoorstroomhoogte"]
            ):
                continue
            elif np.isnan(weir["laagstedoorstroomhoogte"]) or np.isnan(
                weir["hoogstedoorstroomhoogte"]
            ):
                h_laag = h_hoog = np.nanmean(
                    weir[["hoogstedoorstroomhoogte", "laagstedoorstroomhoogte"]].values[:]
                )
            else:
                h_laag = weir["laagstedoorstroomhoogte"]
                h_hoog = weir["hoogstedoorstroomhoogte"]

            geom = weir["geometry"]
            if (not isinstance(geom, Point)) and (branch_geom is not None):
                centroid = geom.centroid
                bool_intersect = branch_geom.intersects(geom)
                if np.sum(bool_intersect.values[:]) > 0:
                    branch = branches_gdf.loc[bool_intersect, "geometry"].iloc[0]
                    branch_loc = branch.project(centroid, normalized=True)
                    branch_loc = np.amax([0.01, np.amin([0.99, branch_loc])])
                    geom = branch.interpolate(branch_loc, normalized=True)
                else:
                    geom = centroid

            if isinstance(weir["soortstuw"], str):
                try:
                    soort_stuw = WEIR_MAPPING[weir["soortstuw"].lower()]
                except KeyError as e:
                    print(e)
                    print("choosing overlaat")
                    soort_stuw = 11

            else:
                soort_stuw = weir["soortstuw"]

            _weir = dict(
                [
                    ("afvoercoefficient", weir["afvoercoefficient_stuw"]),
                    ("code", weir["code"]),
                    ("geometry", geom),
                    ("globalid", weir["globalid"]),
                    ("soortstuw", soort_stuw),
                ]
            )
            weir_list.append(_weir)

            if isinstance(weir["vormopening"], int) | isinstance(weir["vormopening"], float):
                shape = weir["vormopening"]

            elif isinstance(weir["vormopening"], str):  # dtype not a number, so assuming string
                if str(weir["vormopening"]).lower() == "rechthoekig":
                    shape = 3
                else:
                    shape = 1
            else:
                print("Assumed vormopening = 3 for {}".format(ix))
                shape = 3

            opening_gid = str(uuid.uuid4())
            opening = dict(
                [
                    ("afvoercoefficient", weir["afvoercoefficient_opening"]),
                    ("geometry", None),
                    ("globalid", opening_gid),
                    ("hoogstedoorstroombreedte", b_hoog),
                    ("hoogstedoorstroomhoogte", h_hoog),
                    ("laagstedoorstroombreedte", b_laag),
                    ("laagstedoorstroomhoogte", h_laag),
                    ("stuwid", weir["globalid"]),
                    ("vormopening", shape),
                ]
            )
            opening_list.append(opening)

            if isinstance(weir["soortregelbaarheid"], int) | isinstance(
                weir["soortregelbaarheid"], float
            ):
                management = weir["soortregelbaarheid"]

            elif isinstance(
                weir["soortregelbaarheid"], str
            ):  # dtype not a number, so assuming string
                try:
                    management = MANAGEMENT_DEVICE_TYPES[weir["soortregelbaarheid"]]
                except KeyError as e:
                    print(e)
                    print("choosing niet regelbaar")
                    management = 1

            management_device = dict(
                [
                    ("code", weir["code"]),
                    ("geometry", geom),
                    ("globalid", weir["globalid"]),
                    ("kunstwerkopeningid", opening_gid),
                    ("overlaatonderlaat", weir["overlaatonderlaat"]),
                    ("soortregelbaarheid", management),
                    ("stuwid", weir["globalid"]),
                ]
            )
            management_device_list.append(management_device)
        _weir_gdf = gpd.GeoDataFrame(weir_list, geometry="geometry", crs=28992)
        opening_gdf = gpd.GeoDataFrame(opening_list, geometry="geometry", crs=28992)
        management_device_gdf = gpd.GeoDataFrame(
            management_device_list, geometry="geometry", crs=28992
        )

        return _weir_gdf, opening_gdf, management_device_gdf

    def fill_branch_norm_parm_profiles_data(
        defaults, in_branches_gdf: gpd.GeoDataFrame, data_config, insteek_marge=0.25
    ) -> gpd.GeoDataFrame:
        """ """

        def find_branch_width(branch_geom, buffer_list, name, watervlak_geometry):
            for jx, buffer in enumerate(buffer_list):
                buffered_branches = branch_geom.buffer(buffer, cap_style=2)
                polygon = buffered_branches.intersection(watervlak_geometry)

                if jx > 0:
                    overlap_area = polygon.area / buffered_branches.area
                    if overlap_area < 0.9:
                        width = buffer_list[jx - 1] * 2
                        break
                    else:
                        width = buffer * 2

            return dict([("id", name), ("overlap", overlap_area), ("width", width)])

        in_branches_gdf[
            [
                "bodembreedte",
                "bodemhoogte benedenstrooms",
                "bodemhoogte bovenstrooms",
                "hoogte insteek linkerzijde",
                "hoogte insteek rechterzijde",
                "taludhelling linkerzijde",
                "taludhelling rechterzijde",
            ]
        ] = in_branches_gdf[
            [
                "bodembreedte",
                "bodemhoogte benedenstrooms",
                "bodemhoogte bovenstrooms",
                "hoogte insteek linkerzijde",
                "hoogte insteek rechterzijde",
                "taludhelling linkerzijde",
                "taludhelling rechterzijde",
            ]
        ].astype(
            float
        )

        if (
            hasattr(data_config, "peil_gebieden_path")
            and ("diepte" in data_config.np_index_mapping)
            and (data_config.np_index_mapping["diepte"] is not None)
        ):
            peil_gebieden_gdf = gpd.read_file(data_config.peil_gebieden_path)
            peil_gebieden_gdf, _ = map_columns(
                defaults=defaults,
                gdf=peil_gebieden_gdf,
                index_mapping=data_config.peil_index_mapping,
            )

            out_branches_gdf = copy(in_branches_gdf)

            in_branches_with_peil = gpd.sjoin(
                in_branches_gdf, peil_gebieden_gdf, how="left", predicate="intersects"
            )

            for name, branch in in_branches_with_peil.iterrows():
                gid = branch["globalid"]
                bool_ix = out_branches_gdf["globalid"] == gid

                if in_branches_with_peil[in_branches_with_peil["globalid"] == gid].shape[0] > 1:
                    # check for duplicates
                    pass

                if (not check_is_not_na_number(branch["hoogte insteek linkerzijde"])) and (
                    check_is_not_na_number(branch["hoogte insteek rechterzijde"])
                ):
                    out_branches_gdf.loc[
                        bool_ix, "hoogte insteek linkerzijde"
                    ] = out_branches_gdf.at[bool_ix, "hoogte insteek rechterzijde"]
                elif (check_is_not_na_number(branch["hoogte insteek linkerzijde"])) and (
                    not check_is_not_na_number(branch["hoogte insteek rechterzijde"])
                ):
                    out_branches_gdf.loc[
                        bool_ix, "hoogte insteek rechterzijde"
                    ] = out_branches_gdf.at[bool_ix, "hoogte insteek linkerzijde"]
                elif (not check_is_not_na_number(branch["hoogte insteek linkerzijde"])) and (
                    not check_is_not_na_number(branch["hoogte insteek rechterzijde"])
                ):
                    # if check_is_not_na_number(branch["vast peil"]):
                    #     peil = branch["vast peil"]
                    # else:
                    #     peil = np.nanmean([branch["boven peil"], branch["onder peil"]])
                    peil = np.nanmean(branch[["boven peil", "onder peil", "vast peil"]])

                    insteek_hoogte = peil + insteek_marge

                    out_branches_gdf.loc[bool_ix, "hoogte insteek linkerzijde"] = insteek_hoogte
                    out_branches_gdf.loc[bool_ix, "hoogte insteek rechterzijde"] = insteek_hoogte

                if (not check_is_not_na_number(branch["bodemhoogte benedenstrooms"])) and (
                    check_is_not_na_number(branch["bodemhoogte bovenstrooms"])
                ):
                    out_branches_gdf.loc[
                        bool_ix, "bodemhoogte benedenstrooms"
                    ] = out_branches_gdf.at[bool_ix, "bodemhoogte bovenstrooms"]
                elif (check_is_not_na_number(branch["bodemhoogte benedenstrooms"])) and (
                    not check_is_not_na_number(branch["bodemhoogte bovenstrooms"])
                ):
                    out_branches_gdf.loc[
                        bool_ix, "bodemhoogte bovenstrooms"
                    ] = out_branches_gdf.at[bool_ix, "bodemhoogte benedenstrooms"]
                elif (not check_is_not_na_number(branch["bodemhoogte benedenstrooms"])) and (
                    not check_is_not_na_number(branch["bodemhoogte bovenstrooms"])
                ):
                    diepte = branch["diepte"]
                    insteek_hoogte = np.nanmean(
                        out_branches_gdf.loc[
                            bool_ix, ["hoogte insteek linkerzijde", "hoogte insteek rechterzijde"]
                        ].values[:]
                    )

                    if diepte == -99:
                        diepte = np.nan

                    bodem_hoogte = insteek_hoogte - diepte
                    out_branches_gdf.loc[bool_ix, "bodemhoogte benedenstrooms"] = bodem_hoogte
                    out_branches_gdf.loc[bool_ix, "bodemhoogte bovenstrooms"] = bodem_hoogte

            in_branches_gdf = copy(out_branches_gdf)

        if hasattr(data_config, "watervlak_path"):
            watervlak_gdf = gpd.read_file(data_config.watervlak_path, geometry="geometry")
            watervlak_geometry = watervlak_gdf.dissolve(by=None).geometry.values[0]

            buffer_list = [1.25, 2.5, 5, 10, 15, 20, 25, 50, 100]
            # cpus = 8
            # pool = Pool(processes=cpus)
            out_branches_gdf = copy(in_branches_gdf)
            results = []
            save = False
            for ix, (name, branch) in enumerate(in_branches_gdf.iterrows()):
                if ((branch["bodembreedte"] is None) or np.isnan(branch["bodembreedte"])) and (
                    (branch["water_width_index"] is None) or np.isnan(branch["water_width_index"])
                ):
                    if not save:
                        save = True
                    branch_geom = branch.geometry

                    kwds = dict(
                        [
                            ("branch_geom", branch_geom),
                            ("buffer_list", buffer_list),
                            ("name", name),
                            ("watervlak_geometry", watervlak_geometry),
                        ]
                    )
                    # results.append(pool.apply_async(find_brapeilgebiedennch_width, kwds=kwds))
                    results.append(find_branch_width(**kwds))

            for res in results:
                # r.wait()
                # res = r.get()
                out_branches_gdf.at[res["id"], "bodembreedte"] = res["width"]

            if save:
                out_branches_gdf.to_file(data_config.branches_path.replace(".shp", "_ww.shp"))

        if "out_branches_gdf" in locals():
            return out_branches_gdf
        else:
            return in_branches_gdf

    ## START OF FUNCTION
    # ddm = DHydamoDataModel()
    defaults = importlib.import_module("dataset_configs." + defaults)
    data_config = importlib.import_module("dataset_configs." + config)
    code_padding = config[:4] + r"_"  # add prefix of length 4 to all objects with codes
    if code_padding == "rijn":
        code_padding = ""

    gdf_dict = {}
    for key, value in data_config.__dict__.items():
        if not isinstance(value, FMGisFile):
            continue

        gdf = load_geo_file(file_path=value.path, layer=value.name)

        if hasattr(value, "column_selection") and (value.column_selection is not None):
            if isinstance(value.column_selection, dict):
                for s_key, s_value in value.column_selection.items():
                    if isinstance(s_value, list):
                        gdf = gdf.loc[gdf[s_key].isin(s_value), :]
                    else:
                        gdf = gdf.loc[gdf[s_key] == s_value, :]
            else:
                raise NameError("not implemented")

        if (defaults is not None) and hasattr(defaults, value.name):
            default = getattr(defaults, value.name)
        else:
            print("no defaults found for key: {}".format(key))
            default = None
        gdf, _ = map_columns(
            code_pad=code_padding + value.name[:2] + r"_",
            defaults=default,
            gdf=gdf,
            index_mapping=value.column_mapping,
        )
        gdf_dict[value.name] = gdf

    for key, gdf in gdf_dict.items():
        if key == "brug":
            dm = merge_to_dm(dm=dm, feature=key, feature_gdf=create_bridge_data(bridge_gdf=gdf))

        elif key == "duiker" or (key == "afsluitmiddel"): #afsluitmiddel als duiker
            dm = merge_to_dm(dm=dm, feature=key, feature_gdf=create_culvert_data(culvert_gdf=gdf))

        elif key == "gemaal":
            peil_gdf=gpd.read_file(data_config.peilgebieden.path)
            gemaal, pomp, sturing = create_pump_data(pump_gdf=gdf,peil_gdf=peil_gdf)
            features = ["gemaal", "pomp", "sturing"]
            gdfs = [gemaal, pomp, sturing]
            dm = merge_multiple_to_dm(dm=dm, features=features, gdfs=gdfs)

        elif key == "gemeten_profiel":
            (
                profielgroep,
                profiellijn,
                profielpunt,
                ruwheidsprofiel,
            ) = create_measured_profile_data(profile_points_gdf=gdf)
            features = ["profielgroep", "profiellijn", "profielpunt", "ruwheidsprofiel"]
            gdfs = [profielgroep, profiellijn, profielpunt, ruwheidsprofiel]
            dm = merge_multiple_to_dm(dm=dm, features=features, gdfs=gdfs)

        elif key == "normprofiel":
            if hasattr(defaults, "peil"):
                default = defaults.peil
            else:
                default = None
            np_gdf = fill_branch_norm_parm_profiles_data(
                defaults=default, in_branches_gdf=gdf, data_config=data_config
            )
            (
                profielgroep,
                profiellijn,
                profielpunt,
                ruwheidsprofiel,
            ) = create_mp_from_np(branches_gdf=np_gdf)
            features = ["profielgroep", "profiellijn", "profielpunt", "ruwheidsprofiel"]
            gdfs = [profielgroep, profiellijn, profielpunt, ruwheidsprofiel]
            dm = merge_multiple_to_dm(dm=dm, features=features, gdfs=gdfs)

        elif (key == "keringen") or (key == "overiglijnelement"):
            dx = data_config.Models.FM.two_d.dx
            dy = data_config.Models.FM.two_d.dy
            min_length = np.mean([dx, dy])

            if min_length is not None:
                gdf = gdf.drop(gdf[gdf.geometry.length < min_length].index)

            dm = merge_to_dm(dm=dm, feature=key, feature_gdf=gdf)

        elif key == "peilgebieden":
            dm = merge_to_dm(dm=dm, feature=key, feature_gdf=gdf)

        elif (key == "stuw") or (key == "sluis"): # or (key == "afsluitmiddel"):
            if "waterloop" in gdf_dict:
                branches_gdf = gdf_dict["waterloop"]
            else:
                branches_gdf = None
            stuw, kunstwerkopening, regelmiddel = create_weir_data(
                branches_gdf=branches_gdf, weir_gdf=gdf
            )
            features = ["stuw", "kunstwerkopening", "regelmiddel"]
            gdfs = [stuw, kunstwerkopening, regelmiddel]
            dm = merge_multiple_to_dm(dm=dm, features=features, gdfs=gdfs)

        elif key == "waterloop":
            branches_gdf = copy(gdf)
            for ix, branch in gdf.iterrows():
                type_ruwheid = check_roughness(branch)
                branches_gdf.at[ix, "typeruwheid"] = type_ruwheid

            if "peilgebieden" in gdf_dict:
                branches_gdf = add_peil_to_branch(
                    branches_gdf=branches_gdf, peil_gdf=gdf_dict["peilgebieden"]
                )

            if (
                hasattr(data_config, "Peil")
                and hasattr(data_config.Peil, "default_peil")
                and (data_config.Peil.default_peil is not None)
            ):
                branches_gdf = add_default_peil_to_branch(
                    branches_gdf=branches_gdf, default_peil=data_config.Peil.default_peil
                )

            if "bodembreedte" in branches_gdf.columns:
                (
                    waterloop,
                    hydroobject_normgp,
                    normgeparamprofielwaarde,
                ) = create_norm_parm_profiles_v2(branches_gdf=branches_gdf)
                features = ["waterloop", "hydroobject_normgp", "normgeparamprofielwaarde"]
                gdfs = [waterloop, hydroobject_normgp, normgeparamprofielwaarde]
                dm = merge_multiple_to_dm(dm=dm, features=features, gdfs=gdfs)
            else:
                dm = merge_to_dm(dm=dm, feature=key, feature_gdf=branches_gdf)
        else:
            print("Key: {} not implemented".format(key))

    if hasattr(data_config, "river_profile_path") and (data_config.river_profile_path is not None):
        if hasattr(data_config, "river_roughness_path") and (
            data_config.river_roughness_path is not None
        ):
            rough_df = pd.read_csv(data_config.river_roughness_path)
        else:
            rough_df = None

        riv_prof_df = pd.read_csv(data_config.river_profile_path)
        (rivier_profielen, rivier_profielen_data, rivier_ruwheid,) = create_river_profiles(
            branches_gdf=branches_gdf,
            riv_prof_df=riv_prof_df,
            code_padding=code_padding + "wl_",
            rough_df=rough_df,
        )
        features = ["rivier_profielen", "rivier_profielen_data", "rivier_ruwheid"]
        gdfs = [rivier_profielen, rivier_profielen_data, rivier_ruwheid]
        dm = merge_multiple_to_dm(dm=dm, features=features, gdfs=gdfs)
    return dm

