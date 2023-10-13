import json
import urllib.parse
import urllib.request
from pathlib import Path

import geopandas as gpd
import requests


def wfs_download_rijnland(url: str, output_file_path: str = None) -> None:
    """ """
    gdf = gpd.read_file(url)

    if output_file_path is not None:
        Path(output_file_path).parents[0].mkdir(exist_ok=True)
        gdf.to_file(output_file_path)
    else:
        return gdf


def json_download_rijnland(
    url: str, output_file_path: str = None, geometry_type: str = "esriGeometryPolygon"
) -> None:
    """ """

    params = {
        "where": "1=1",
        "geometryType": geometry_type,
        "spatialRel": "esriSpatialRelIntersects",
        "relationParam": "",
        "outFields": "*",
        "returnGeometry": "true",
        "geometryPrecision": "",
        "outSR": "",
        "returnIdsOnly": "false",
        "returnCountOnly": "false",
        "orderByFields": "",
        "groupByFieldsForStatistics": "",
        "returnZ": "false",
        "returnM": "false",
        "returnDistinctValues": "false",
        "f": "pjson",
    }

    encode_params = urllib.parse.urlencode(params).encode("utf-8")

    response = urllib.request.urlopen(url, encode_params)
    json_response = response.read()

    Path(output_file_path).parents[0].mkdir(exist_ok=True)
    json_file_path = Path(output_file_path).parents[0] / "temp.json"

    with open(json_file_path, "wb") as ms_json:
        ms_json.write(json_response)

    # print(json_response)

    gdf = gpd.read_file(json_file_path)

    if output_file_path is not None:
        gdf.to_file(output_file_path)
    else:
        return gdf
