import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import LineString, MultiLineString

input_path = r"D:\Work\Project\HL-23006\GIS\Keringen_met_hoogte\hdsr.shp"
output_path = r"D:\Work\Project\HL-23006\GIS\Keringen_met_hoogte\hdsr_simplified.shp"


def split_line_by_point(line, max_length: float = 500, max_seg_length: float = 25):
    length = line.length
    if length < max_length:
        return line

    n_seg = np.ceil(length / max_length).astype(int)
    l_list = []
    for n in range(0, n_seg):
        sp = n / n_seg
        ep = (n + 1) / n_seg
        points_norm = np.linspace(sp, ep, np.ceil(max_length / max_seg_length).astype(int))
        points = []
        for p in points_norm:
            points.append(line.interpolate(p, normalized=True))
        l_list.append(LineString(points))

    # mp = MultiPoint(p_list)
    mls = MultiLineString(l_list)
    return mls


gdf = gpd.read_file(input_path)
print(gdf.head())




result = (
    gdf.assign(
        geometry=gdf.apply(
            lambda x: split_line_by_point(
                x.geometry,
            ),
            axis=1,
        )
    )
    .explode(index_parts=False)
    .reset_index(drop=True)
)
result["geometry"] = result["geometry"].simplify(tolerance=1)

result["length"] = result["geometry"].length

result.to_file(output_path)
