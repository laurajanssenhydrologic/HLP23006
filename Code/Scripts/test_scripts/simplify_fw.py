import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import LineString, MultiLineString

gdf_1 = gpd.read_file(
    r"D:\Work\Project\P1414\GIS\HYDAMO\Combined_test_v20_WBD.gpkg", layer="keringen"
)
gdf_2 = gpd.read_file(
    r"D:\Work\Project\P1414\GIS\HYDAMO\Combined_test_v20_WBD.gpkg", layer="overiglijnelement"
)


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


gdf = gpd.GeoDataFrame(pd.concat([gdf_1, gdf_2], ignore_index=True), crs=gdf_1.crs)
gdf.to_file(r"D:\Work\Project\P1414\GIS\Keringen\fw_orig.shp")
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

gdf.to_file(r"D:\Work\Project\P1414\GIS\Keringen\fw_test_1.shp")
result.to_file(r"D:\Work\Project\P1414\GIS\Keringen\fw_test_2.shp")
