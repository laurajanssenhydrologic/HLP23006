import uuid
from typing import Tuple

import geopandas as gpd
import numpy as np
from hydrolib.core.io.ext.models import ExtModel, Lateral
from hydrolib.dhydamo.core.drr import DRRModel
from scipy.spatial import KDTree
from shapely.geometry import Point
from tqdm import tqdm


def connect_external_weirs_to_branches(
    knopen_gdf: gpd.GeoDataFrame,
    wl_gdf: gpd.GeoDataFrame,
    drrmodel: DRRModel = None,
    index_colum: str = "CODE",
) -> Tuple[DRRModel, ExtModel]:
    if drrmodel is None:
        drrmodel = DRRModel()

    assert isinstance(drrmodel, DRRModel)

    wl_gdf["gid"] = [str(uuid.uuid4()) for _ in range(wl_gdf.shape[0])]
    wl_point_list = []
    wl_list = []
    for ix, wl in wl_gdf.iterrows():
        coords = wl.geometry.coords[:]
        gid = wl["gid"]
        for x, y in coords:
            wl_point_list.append([x, y])
            wl_list.append(gid)

    wl_len = len(wl_list)
    wl_kdtree = KDTree(data=wl_point_list)

    lateral_list = []
    for ix, knoop in tqdm(knopen_gdf.iterrows(), total=knopen_gdf.shape[0]):
        # Create rr paved node
        b_node_name = str(knoop[index_colum]) + "_lat"
        area = knoop["verhard_m2"]
        if float(area) < 1e-2:
            tqdm.write(b_node_name)
        drrmodel.paved.add_paved(
            id=knoop[index_colum],
            area=area,
            surface_level=knoop["ahn_mean"],
            street_storage=10,
            sewer_storage=knoop["berging_mm"],
            pump_capacity=knoop["geinst_cap"],
            meteo_area=knoop[index_colum],
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
        wl_geom = wl_gdf.loc[wl_gdf["gid"] == wl_id, "geometry"].values[0]
        wl_name = wl_gdf.loc[wl_gdf["gid"] == wl_id, "CODE"].values[0]
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

        drrmodel.external_forcings.add_boundary_node(id=b_node_name, px=px, py=py)

    extforcefilenew = ExtModel(lateral=lateral_list)
    return drrmodel, extforcefilenew


if __name__ == "__main__":
    knopen_path = r"D:\Work\Project\HL-23006\GIS\Riool\comb.shp"
    wl_path = r"D:\Work\Project\HL-23006\GIS\Selectie_watergangen_21062023\sanitized.shp"

    knopen_gdf = gpd.read_file(knopen_path)
    knopen_gdf["berging_mm"] = knopen_gdf["berging_mm"].fillna(value=0)
    knopen_gdf["geinst_cap"] = knopen_gdf["geinst_cap"].fillna(value=0)
    knopen_gdf["verhard_m2"] = knopen_gdf["verhard_m2"].fillna(value=0)

    wl_gdf = gpd.read_file(wl_path)

    drrmodel, extforcefilenew = connect_external_weirs_to_branches(
        drrmodel=None, index_colum="objectid_1", knopen_gdf=knopen_gdf, wl_gdf=wl_gdf
    )
