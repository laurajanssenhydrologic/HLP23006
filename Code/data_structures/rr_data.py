import importlib
import sys
import uuid
from copy import deepcopy
from datetime import datetime, timedelta
from pathlib import Path
from typing import List, Tuple

import geopandas as gpd
import numpy as np
import xarray as xr
from hydrolib.core.io.bc.models import ForcingModel, QuantityUnitPair, TimeSeries
from hydrolib.core.io.dimr.models import RRComponent
from hydrolib.core.io.ext.models import ExtModel, Lateral
from hydrolib.dhydamo.core.drr import DRRModel
from hydrolib.dhydamo.io.drrwriter import DRRWriter
from scipy.spatial import KDTree
from shapely.geometry import Point
from tqdm import tqdm

from data_structures.fm_data import FMData


class RRData:
    def __init__(self):
        self.rr: RRComponent = None
        self.drrmodel: DRRModel = None
        self.fm: FMData = None

    def connect_external_weirs_to_branches(
        self,
        knopen_gdf: gpd.GeoDataFrame,
        wl_gdf: gpd.GeoDataFrame,
        knoop_index_column: str = "CODE",
        wl_index_column: str = "CODE",
    ) -> Tuple[DRRModel, ExtModel]:
        if self.drrmodel is None:
            self.drrmodel = DRRModel()

        assert isinstance(self.drrmodel, DRRModel)

        # wl_gdf["gid"] = [str(uuid.uuid4()) for _ in range(wl_gdf.shape[0])]
        wl_point_list = []
        wl_list = []
        for ix, wl in wl_gdf.iterrows():
            coords = wl.geometry.coords[:]
            gid = wl["globalid"]
            for x, y in coords:
                wl_point_list.append([x, y])
                wl_list.append(gid)

        wl_len = len(wl_list)
        wl_kdtree = KDTree(data=wl_point_list)

        knopen_gdf["berging_mm"] = knopen_gdf["berging_mm"].fillna(value=0)
        knopen_gdf["geinst_cap"] = knopen_gdf["geinst_cap"].fillna(value=0)
        knopen_gdf["verhard_m2"] = knopen_gdf["verhard_m2"].fillna(value=0)

        lateral_list = []
        for ix, knoop in tqdm(knopen_gdf.iterrows(), total=knopen_gdf.shape[0]):
            # Create rr paved node
            b_node_name = str(knoop[knoop_index_column]) + "_lat"
            area = knoop["verhard_m2"]
            if float(area) < 1e-2:
                continue

            self.drrmodel.paved.add_paved(
                id=knoop[knoop_index_column],
                area=area,
                surface_level=knoop["ahn_mean"],
                street_storage=10,
                sewer_storage=knoop["berging_mm"],
                pump_capacity=knoop["geinst_cap"],
                meteo_area=knoop[knoop_index_column],
                px=knoop.geometry.x,
                py=knoop.geometry.y,
                boundary_node=b_node_name,
            )
            # Find nearest point in branches
            coords = np.array(knoop.geometry.coords[:])[0]
            _, ix_point_in_wl = wl_kdtree.query(coords)
            point_in_wl = wl_point_list[ix_point_in_wl]
            wl_id = wl_list[ix_point_in_wl]

            point = Point(point_in_wl)
            wl_geom = wl_gdf.loc[wl_gdf["globalid"] == wl_id, "geometry"].values[0]
            wl_name = wl_gdf.loc[wl_gdf["globalid"] == wl_id, wl_index_column].values[0]
            chainage = wl_geom.project(point)
            wl_point = wl_geom.interpolate(chainage)
            px = wl_point.x
            py = wl_point.y

            lateral = Lateral(
                id=b_node_name,
                name=b_node_name,
                branchId=wl_name,
                chainage=chainage,
                discharge="realtime",
            )
            lateral_list.append(lateral)

            self.drrmodel.external_forcings.add_boundary_node(id=b_node_name, px=px, py=py)

        extforcefilenew = ExtModel(lateral=lateral_list)
        self.fm.merge_ext_files(extforcefile=extforcefilenew)
        # return extforcefilenew

    def convert_MF_MS_to_lat(
        self,
        ae_gdf: gpd.GeoDataFrame,
        laterals: List,
        network: xr.Dataset,
        d_ref: datetime = datetime(2001, 1, 1, 0, 0, 0),
    ) -> Tuple[ExtModel, ForcingModel]:
        print("starting conversion")
        # ae_gdf = gpd.read_file(ae_path)

        # Load mesh nodes in KDtree to find branchid and chainage for
        branches_ix = network["mesh1d_node_branch"].values
        branches_id = network["network_branch_id"].values
        offsets = network["mesh1d_node_offset"].values
        xs = network["mesh1d_node_x"].values
        ys = network["mesh1d_node_y"].values

        coords = np.stack([xs, ys], axis=0).T
        kdtree = KDTree(data=coords)

        _laterals = []
        _forcings = []
        # forcing_file = "boundaryconditions.bc"
        ## Convert to HUDROLIB core
        for lateral in tqdm(laterals):
            name = lateral[0]
            name = name.removeprefix("'drain").removesuffix("'")
            ae = ae_gdf.loc[ae_gdf["CODE_ORIG"] == name, :]

            centre = ae.centroid
            x, y = centre.x.values[0], centre.y.values[0]
            # print([x, y])
            dist, nearest = kdtree.query([x, y], k=1)
            branchix = branches_ix[nearest]
            branchid = branches_id[branchix]
            offset = offsets[nearest]
            # print(name, branchid, offset)

            _forcing = TimeSeries(
                datablock=lateral[1:],
                name=name,
                timeinterpolation="linear",
                quantityunitpair=[
                    QuantityUnitPair(
                        quantity="time",
                        unit="hours since {}".format(d_ref.strftime(r"%Y-%m-%d %H:%M:%S")),
                    ),
                    QuantityUnitPair(quantity="lateral_discharge", unit="m³/s"),
                ],
            )
            _forcings.append(_forcing)
            _lateral = Lateral(
                id=name,
                name=name,
                branchId=branchid,
                chainage=offset,
                discharge=ForcingModel(forcing=_forcing),
            )
            _laterals.append(_lateral)

        # extforcefilenew = ExtModel(lateral=_laterals)
        extforcefilenew = ExtModel()
        extforcefilenew.Config.validate_assignment = False
        extforcefilenew.lateral = _laterals
        forcingmodel = ForcingModel(forcing=_forcings)
        # forcingmodel = ForcingModel()
        # forcingmodel.Config.validate_assignment = False
        # forcingmodel.forcing = _forcings
        self.fm.files_to_save["boundaryconditions.bc"] = forcingmodel

        self.fm.merge_ext_files(extforcefile=extforcefilenew)

        # return extforcefilenew, forcingmodel

    def read_MF_MS_output(
        self,
        lat_path: str,
        date_format: str = r"'%Y/%m/%d;%H:%M:%S'",
        d_ref: datetime = datetime(2001, 1, 1, 0, 0, 0),
        total: int = 14774565,
    ) -> List:
        print("reading file")
        ## Read as text file and parse headers and data including timestamp
        laterals = []
        entries = 0
        with open(lat_path) as file:
            for n, line in tqdm(enumerate(file), total=total):
                if "FLBR" in line:
                    sline = line.split()
                    entry = [sline[2]]
                elif "TBLE" in line:
                    continue
                elif "tble flbr" in line:
                    entries += 1
                    laterals.append(entry)
                    # if entries > 1:
                    #     break
                    # break
                else:
                    date, value, _ = line.split()
                    _date = datetime.strptime(date, date_format) - d_ref
                    date = str(_date.total_seconds() // 3600)
                    line = [date, value]
                    entry.append(line)

        return laterals

    def set_basemodel_from_config(self, config: str):
        model_config = getattr(importlib.import_module("dataset_configs." + config), "Models")
        if not hasattr(model_config, "RR"):
            print("No RR config found")
            return

        if self.drrmodel is None:
            self.drrmodel = DRRModel()

        start_time = model_config.RR.start_time
        stop_time = model_config.RR.stop_time
        timestep = model_config.RR.timestep

        start_time = datetime.strptime(str(start_time), "%Y%m%d")
        stop_time = start_time + timedelta(seconds=stop_time)

        self.drrmodel.d3b_parameters["Timestepsize"] = timestep
        self.drrmodel.d3b_parameters["StartTime"] = "'{:0>4d}/{:0>2d}/{:0>2d};00:00:00'".format(
            start_time.year, start_time.month, start_time.day
        )  # should be equal to refdate for D-HYDRO
        self.drrmodel.d3b_parameters["EndTime"] = "'{:0>4d}/{:0>2d}/{:0>2d};00:00:00'".format(
            stop_time.year, stop_time.month, stop_time.day
        )
        self.drrmodel.d3b_parameters["RestartIn"] = 0
        self.drrmodel.d3b_parameters["RestartOut"] = 0
        self.drrmodel.d3b_parameters["RestartFileNamePrefix"] = "Test"
        self.drrmodel.d3b_parameters["UnsaturatedZone"] = 1
        self.drrmodel.d3b_parameters["UnpavedPercolationLikeSobek213"] = -1
        self.drrmodel.d3b_parameters["VolumeCheckFactorToCF"] = 100000

    def write_drr(self, output_path: Path, wwtp: Tuple = (135888, 458035)):
        if not hasattr(self, "drrmodel"):
            print("No drrmodel model present")
            return None
        else:
            drr = self.drrmodel

        print("Saving RR model to disk")
        rr_writer = DRRWriter(rrmodel=drr, output_dir=output_path, name="RR", wwtp=wwtp)
        rr_writer.write_all()
        return drr
