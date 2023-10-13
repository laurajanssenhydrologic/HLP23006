import importlib
from typing import List, Union

import geopandas as gpd
import numpy as np
import pandas as pd
from geo_tools.clip_tools import _clip_structures_by_branches
from geo_tools.merge_networks import merge_networks

from data_structures.fixedweirs_helpers import (create_dambreak_data,
                                                create_fixed_weir_data)
from data_structures.hydamo_helpers import (check_and_fix_duplicate_code,
                                            convert_to_dhydamo_data)
from data_structures.roi_data_model import FMDataModel as DataModel
from data_structures.to_dhydro_helpers import fm_to_dhydro, write_dimr


class DHydroData:
    """
    Helper class to convert raw data to:
    1. load raw data into a DHydamoDataModel
    2. load a DHydamo geopackage into a DHydamoDataModel
    3. save a DHydamoDataModel into a DHydamo geopackage
    4. save a DHYdamoDataModel into a D-HYDRO model
    5. Clip structures by (buffered) branches extent

    Attributes:
        dm (DHydamoDataModel): instance of DHydamoDataModel that is used by this clas
        features (list): list of features that are present in self.dm and do not have None values
    """

    def __init__(self):
        pass

    def clip_structures_by_branches(
        self, buffer: float = 1, features: list = None, min_overlap: float = 0.5
    ) -> None:
        """
        Class method to drop structures from the data that are farther from a branch than "buffer"
        Structures that are linestring are required to overlap with a branch for at least "min_overlap"

        Args:
            buffer (float): distance from branch that structures are kept (default = 1m)
            min_overlap (float): percentage of branch and LineString structure that should overlap in order to keep it.
                                 this drop, for example, culverts that are perpendicular to a branch

        Returns:
            None
        """
        if (not hasattr(self, "dm")) or (not hasattr(self, "features")):
            raise AttributeError("Database not loaded")

        if features is None:
            features = self.features

        self.dm = _clip_structures_by_branches(
            dm=self.dm, features=features, buffer=buffer, min_overlap=min_overlap
        )

    def check_missing_values(
        self, missing_values: List = [None, np.nan, -999, 999, -99, 99, "n.v.t."], fill_value=None
    ):
        dm = DataModel()
        dm.Config.validate_assignment = False

        invalid_objs = []
        for key, value in self.dm.__dict__.items():
            # check if attribute is not None
            if value is None:
                continue

            gdf = value
            _gdf = gdf.drop(columns="geometry").replace(
                to_replace=missing_values, value=fill_value
            )
            invalid_rows = gdf[_gdf.isnull().any(axis=1)].replace(
                to_replace=missing_values, value=fill_value
            )
            setattr(dm, key, invalid_rows)

        return dm

    def dambreaks_from_config(self, config: str, defaults: str):
        dm = DataModel()
        dm = create_dambreak_data(config=config, defaults=defaults, dm=dm)

        self._set_dm(dm=dm)

    def delete_incomplete_features(
        self,
        features: List,
        columns: List[List],
        missing_values: List = [None, np.nan, -999, 999, -99, 99, "n.v.t."],
    ):
        if (not hasattr(self, "dm")) or (not hasattr(self, "features")):
            raise AttributeError("Database not loaded")

        if len(features) != len(columns):
            raise ValueError("both lists should be the same size")

        for feature, _columns in zip(features, columns):
            if not feature in self.features:
                raise AttributeError("Unknown feature(s) " + feature)

            gdf = getattr(self.dm, feature)
            nentries = gdf.shape[0]
            if isinstance(_columns, list):
                for column in _columns:
                    gdf = gdf.loc[~gdf[column].isin(missing_values), :]
            else:
                gdf = gdf.loc[~gdf[_columns].isin(missing_values), :]

            print("dropped {} entries from {}".format(nentries - gdf.shape[0], feature))
            setattr(self.dm, feature, gdf)

    def fixed_weirs_from_raw_data(self, config: str, defaults: str, min_length: float = None):
        dm = DataModel()
        dm = create_fixed_weir_data(config=config, defaults=defaults, dm=dm, min_length=min_length)
        self._set_dm(dm=dm)

    def hydamo_from_raw_data(
        self,
        config: str,
        defaults: str,
        branch_snap_dist: float = 10,
        validate: bool = True,
    ) -> None:
        """
        Class method to load raw_data into a DHydamoDataModel. This datamodel validates data against expected values

        Args:
            defaults (str): default settings to use (should be in ./dataset_configs)
            config (str): configuration file to use (should be in ./dataset_configs)

        Returns:
            None
        """
        # load features and add to DHydamoDataModel
        dm = DataModel()
        dm.Config.validate_assignment = validate
        dm = convert_to_dhydamo_data(dm=dm, defaults=defaults, config=config)

        self._set_dm(branch_snap_dist=branch_snap_dist, dm=dm)

    def hydamo_from_gpkg(
        self,
        gpkg_path: str,
        branch_snap_dist: float = 10,
    ) -> None:
        """
        Class method to load data from geopackage and validate data against DHydamoDataModel

        Args:
            gpkg_path (str): file location of the geopackage to load

        Returns
            None
        """

        # set geopackage path to self
        # self.gpkg_path = gpkg_path

        # initialize datamodel
        dm = DataModel()
        # loop over datamodel attributes and check if they are presentin the geopackage
        # if so, set them to the datamodel
        attributes = dm.__dict__.keys()
        for attribute in attributes:
            try:
                # print("succesfully loaded {}".format(attribute))
                data = gpd.read_file(gpkg_path, layer=attribute)
            except ValueError:
                # print("failed to load {}".format(attribute))
                continue

            setattr(dm, attribute, data)

        # set datamodel to self
        self._set_dm(branch_snap_dist=branch_snap_dist, dm=dm)

    def hydamo_to_gpkg(self, output_gpkg: str) -> None:
        """
        Class method that saves DHydamoDataModel to geopackage

        Args:
            output_gpkg (str): file location to save geopackage to

        Returns:
            None
        """
        # self.gpkg_path = output_gpkg
        self.dm.to_gpkg(output_gpkg=output_gpkg)

    def to_dhydro(self, config: str, output_folder: str, defaults: str = "defaults", write=True):
        """
        Class method that converts a DHydamoDataModel to a D-HYDRO Model and saves unless write=False

        Args:
            config (str): configuration file to use (should be in ./dataset_configs)
            output_folder (str): folder to save D-HYDRO model to

        Returns:
            None
        """
        # check if data has been loaded and correct attributes are set
        if (not hasattr(self, "dm")) | (not hasattr(self, "features")):
            raise AttributeError("Modeldatabase not loaded")

        # if hasattr(importlib.import_module("dataset_configs." + config), "Dambreak"):
        #     self.dambreaks_from_config(config=config, defaults=defaults)
        fm_to_dhydro(self=self, config=config, output_folder=output_folder)

        if write:
            self.write_dimr(output_folder=output_folder)

    def validate_dm(self):
        dm = DataModel()
        dm.Config.validate_assignment = True

        for key, value in self.dm.__dict__.items():
            # check if attribute is not None
            if value is not None:
                try:
                    setattr(dm, key, value)
                except Exception as e:
                    print(e)

    def write_dimr(self, output_folder: str):
        return write_dimr(fm=self.fm, output_folder=output_folder)

    def _set_dm(self, dm: DataModel, branch_snap_dist: float = None) -> None:
        """
        Class method to add a DHydamoDataModel to self while checking for pre-existing data.
        This method simply assigns the data if there is no pre-existing data.
        Otherwise, it appends to the pre-existing data

        Args:
            dm (DHydamoDataModel): a DHydamoDataModel that is to be added to self.dm

        Returns:
            None
        """

        # check if a datamodel already exists
        if not hasattr(self, "dm"):
            # if not, just assign
            self.dm = dm

            # and fill list of features present
            self.features = []
            for key, value in self.dm.__dict__.items():
                # check if attribute is not None
                if value is not None:
                    self.features.append(key)

        # if it does exist, check per attribute in datamodel if it already exists
        else:
            for key, value in dm.__dict__.items():
                # check if attribute is not None
                if value is not None:
                    # if an attribute does not exist, simply assign and update features list
                    if getattr(self.dm, key) is None:
                        setattr(self.dm, key, getattr(dm, key))
                        self.features.append(key)

                    # else, concatonate the old and new geodataframe and assign to the datamodel
                    else:
                        in_gdf = getattr(dm, key)

                        if key == "waterloop":
                            in_gdf = merge_networks(
                                data_base_input=in_gdf,
                                data_match_input=getattr(self.dm, key),
                                max_dist=branch_snap_dist,
                            )

                        new_gdf = gpd.GeoDataFrame(
                            pd.concat([getattr(self.dm, key), in_gdf], ignore_index=True)
                        )
                        new_gdf = check_and_fix_duplicate_code(new_gdf)

                        setattr(self.dm, key, new_gdf)
