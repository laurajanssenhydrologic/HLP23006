from pathlib import Path
from typing import Union

import geopandas as gpd
import xugrid as xu
from hydrolib.core.io.mdu.models import FMModel
from shapely.geometry import Point, Polygon
from tqdm import tqdm


def bounds_from_grid(mapfile: str, output_file: str = None) -> gpd.GeoDataFrame:
    ds = xu.open_dataset(mapfile)
    mesh2d = ds.ugrid.grids[2]

    fnc = mesh2d.face_node_connectivity
    nx = mesh2d.node_x
    ny = mesh2d.node_y

    cell_list = []
    for ix in tqdm(range(fnc.shape[0])):
        nodes = fnc[ix, :]
        point_list = []
        for node in nodes:
            point = Point(nx[node], ny[node])
            point_list.append(point)
        polygon = Polygon(point_list)
        cell = dict([("ix", ix), ("geometry", polygon)])
        cell_list.append(cell)

    gdf = gpd.GeoDataFrame(cell_list).dissolve(by=None)

    if output_file is not None:
        gdf.to_file(output_path)

    return gdf


def _set_modelrun_files(mdupath: Union[Path, str]):
    """Derive output file paths for a given D-Flow FM mdu file.

    Args:
        mdupath (Union[Path, str]): Path to the D-Flow FM MDU input file.
    """
    # global modelrun_mdufile, modelrun_diafile, modelrun_mapfile, modelrun_clmfile, modelrun_hisfile, fmmodel

    if isinstance(mdupath, str):
        modelrun_mdufile = Path(mdupath)
    elif isinstance(mdupath, Path):
        modelrun_mdufile = mdupath
    else:
        raise ValueError("Argument mdupath must be either str or Path.")

    md_ident = modelrun_mdufile.stem
    workdir = modelrun_mdufile.parent

    fmmodel = FMModel(modelrun_mdufile, recurse=False)
    outputdir = Path(fmmodel.output.outputdir or "DFM_OUTPUT_" + md_ident)
    if not outputdir.is_absolute():
        outputdir = workdir / outputdir

    modelrun_diafile = outputdir / (md_ident + ".dia")
    modelrun_mapfile = outputdir / (fmmodel.output.mapfile.filepath or md_ident + "_map.nc")
    modelrun_clmfile = outputdir / (fmmodel.output.classmapfile.filepath or md_ident + "_clm.nc")
    modelrun_hisfile = outputdir / (fmmodel.output.hisfile.filepath or md_ident + "_his.nc")

    print(f"Using model files:")
    print(f" - MDU file : {modelrun_mdufile}")
    print(f" - .dia file: {modelrun_diafile}")
    print(f" - map file : {modelrun_mapfile}")
    print(f" - clm file : {modelrun_clmfile}")
    print(f" - his file : {modelrun_hisfile}")

    return modelrun_diafile, modelrun_clmfile, modelrun_hisfile, modelrun_mapfile


if __name__ == "__main__":
    mapfile = r"D:\Work\Project\P1414\Models\Combined\V21_WBD_HD\run_model_changed_coords\dflowfm\output\DFM_map.nc"
    output_path = r"D:\Work\Project\P1414\Models\Combined\V21_WBD_HD\run_model_changed_coords\dflowfm\output\clip.shp"
    gdf = bounds_from_grid(mapfile=mapfile, output_file=output_path)
