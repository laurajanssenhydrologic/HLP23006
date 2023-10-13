from datetime import datetime
from typing import List, Tuple

import geopandas as gpd
import numpy as np
import xarray as xr
from hydrolib.core.io.bc.models import ForcingModel, QuantityUnitPair, TimeSeries
from hydrolib.core.io.ext.models import ExtModel, Lateral
from scipy.spatial import KDTree
from tqdm import tqdm


def read_MF_MS_output(
    lat_path: str,
    date_format: str = r"'%Y/%m/%d;%H:%M:%S'",
    d_ref: datetime = datetime(2001, 1, 1, 0, 0, 0),
) -> List:
    laterals = []
    entries = 0
    with open(lat_path) as file:
        for n, line in tqdm(enumerate(file), total=14774565):
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


def convert_MF_MS_to_lat(
    ae_gdf: gpd.GeoDataFrame,
    laterals: List,
    network: xr.Dataset,
    d_ref: datetime = datetime(2001, 1, 1, 0, 0, 0),
) -> Tuple[ExtModel, ForcingModel]:
    print("reading file")
    ## Read as text file and parse headers and data including timestamp

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

    extforcefilenew = ExtModel(lateral=_laterals)
    forcingmodel = ForcingModel(forcing=_forcings)

    return extforcefilenew, forcingmodel


if __name__ == "__main__":
    ae_path = r"D:\Work\Project\HL-23006\GIS\koppeling ms sobek\Afwateringseenheden_REFERENTIE_2021-11-12.geojson"
    net_path = r"D:\Work\Project\HL-23006\models\HDSR\V0\export\dflowfm\network.nc"
    ext_path = r"D:\Work\Project\HL-23006\GIS\koppeling ms sobek\test_v2.ext"
    forcing_path = r"D:\Work\Project\HL-23006\GIS\koppeling ms sobek\boundaryconditions_v2.bc"
    lat_path = r"D:\Work\Project\HL-23006\GIS\koppeling ms sobek\LATERAL.DAT"

    ae_gdf = gpd.read_file(ae_path)
    laterals = read_MF_MS_output(lat_path=lat_path)
    network = xr.open_dataset(net_path)

    extforcefilenew, forcingmodel = convert_MF_MS_to_lat(
        ae_gdf=ae_gdf, laterals=laterals, network=network
    )

    extforcefilenew.save(filepath=ext_path, recurse=True)
    forcingmodel.save(filepath=forcing_path)
