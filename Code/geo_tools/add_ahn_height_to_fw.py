from copy import copy

import geopandas as gpd
import networkx as nx
import numpy as np
import rasterio
from rasterstats import zonal_stats
from shapely.geometry import LineString, Point


def add_height_to_linestrings(
    gdf: gpd.GeoDataFrame, ahn_path: str, buffer: float = None, n_std: float = 2
) -> gpd.GeoDataFrame:
    ahn = rasterio.open(ahn_path)
    affine = ahn.transform
    array = ahn.read(1)
    out_gdf = copy(gdf)
    if buffer is None:
        for name, row in gdf.iterrows():
            line = row.geometry
            _line = np.asarray(line.coords)[:, 0:2]
            # np_line = np.insert(_line, 2, values=np.nan, axis=1)

            heights = ahn.sample(_line)

            z = np.array([z for z in heights]).flatten()
            z_mean = np.nanmean(z)
            z[np.isnan(z)] = z_mean

            if np.isnan(z_mean):
                out_gdf.drop(name, inplace=True)
                print("dropped {}".format(name))
                continue

            np_line = np.insert(_line, 2, values=z, axis=1)
            # np_line[:, 2] = [h for h in heights]

            # h_mean = np.nanmean(np_line[:, 2])
            # np_line[np.isnan(np_line[:, 2]), 2] = h_mean
            out_gdf.loc[name, "z_min"] = np.nanmin(z)
            out_gdf.loc[name, "geometry"] = LineString(np_line)

    else:

        for name, row in gdf.iterrows():
            print(name)
            line = row.geometry
            _line = np.asarray(line.coords)[:, 0:2]
            points = []
            for x, y in _line:
                point = Point(x, y)
                b_point = point.buffer(buffer)
                points.append(b_point)
            # obtain mean height around all points of a line
            stats = zonal_stats(
                points,
                array,
                affine=affine,
                stats="",
                add_stats={"nanmax": np.nanmax},
                nodata=np.nan,
            )
            heights = []
            for result in stats:
                heights.append(result["nanmax"])

            # convert to numpy array and remove values more than two standard devations from the nanmean
            z = np.array([z for z in heights]).flatten()
            z_mean = np.nanmean(z)
            z_std = np.nanstd(z)
            bool_arry = (z < (z_mean - n_std * z_std)) | (z > (z_mean + n_std * z_std))
            if np.sum(bool_arry) == z.shape[0]:
                z[:] = z_mean
            else:
                z[bool_arry] = np.nan
                # Fill all nan-values with the new median
                z[np.isnan(z)] = np.nanmedian(z)

            if np.isnan(z_mean):
                out_gdf.drop(name, inplace=True)
                print("dropped {}".format(name))
                continue

            np_line = np.insert(_line, 2, values=z, axis=1)
            # np_line[:, 2] = [h for h in heights]

            # h_mean = np.nanmean(np_line[:, 2])
            # np_line[np.isnan(np_line[:, 2]), 2] = h_mean
            out_gdf.loc[name, "z_min"] = np.nanmin(z)
            out_gdf.loc[name, "geometry"] = LineString(np_line)
    ahn.close()
    return out_gdf


if __name__ == "__main__":
    ahn_path = r"D:\Work\Project\P1414\GIS\AHN\AHN4_WSS_filled.tif"
    in_out_dict = dict(
        [
            (
                r"D:\Work\Project\P1414\GIS\Keringen\Noordzeekeringen.shp",
                r"D:\Work\Project\P1414\GIS\Keringen_met_hoogte\Noordzeekeringen.shp",
            ),
            (
                r"D:\Work\Project\P1414\GIS\RWS\oeverconstructie_verticaal_clipped.shp",
                r"D:\Work\Project\P1414\GIS\Keringen_met_hoogte\ARK.shp",
            ),
            (
                "D:\Work\Project\P1414\GIS\HDSR\Legger\Waterkeringen_Lijn\Waterkeringen_Lijn_BR.shp",
                "D:\Work\Project\P1414\GIS\Keringen_met_hoogte\hdsr.shp",
            ),
            (
                r"D:\Work\Project\P1414\GIS\HHDelfland\Legger_Delfland_shp\Waterkeringen\Zeewering landwaartse begrenzing.shp",
                r"D:\Work\Project\P1414\GIS\Keringen_met_hoogte\hhd_zeewering.shp",
            ),
            (
                r"D:\Work\Project\P1414\GIS\HHDelfland\Legger_Delfland_shp\Waterkeringen\Regionale waterkering buitenkruinlijn.shp",
                r"D:\Work\Project\P1414\GIS\Keringen_met_hoogte\hhd_regionale_kering.shp",
            ),
            (
                r"D:\Work\Project\P1414\GIS\HHDelfland\Legger_Delfland_shp\Waterkeringen\Overige waterkering polderkade middenkruinlijn.shp",
                r"D:\Work\Project\P1414\GIS\Keringen_met_hoogte\hhd_polderkade.shp",
            ),
            (
                r"D:\Work\Project\P1414\GIS\HHDelfland\Legger_Delfland_shp\Waterkeringen\Overige waterkering landscheiding middenkruinlijn.shp",
                r"D:\Work\Project\P1414\GIS\Keringen_met_hoogte\hhd_landscheiding.shp",
            ),
            (
                r"D:\Work\Project\P1414\GIS\HHDelfland\Legger_Delfland_shp\Waterkeringen\Delflandsedijk buitenkruinlijn.shp",
                r"D:\Work\Project\P1414\GIS\Keringen_met_hoogte\hhd_delflandsedijk.shp",
            ),
            (
                r"D:\Work\Project\P1414\GIS\HHRijnland\Primaire_keringen\Primaire_kering.shp",
                r"D:\Work\Project\P1414\GIS\Keringen_met_hoogte\hhr_primaire_kering.shp",
            ),
            (
                r"D:\Work\Project\P1414\GIS\HHRijnland\Regionale_keringen\Regionale_keringen.shp",
                r"D:\Work\Project\P1414\GIS\Keringen_met_hoogte\hhr_regionale_kering.shp",
            ),
            (
                r"D:\Work\Project\P1414\GIS\HHSK\Keringen\Primaire_waterkering.shp",
                r"D:\Work\Project\P1414\GIS\Keringen_met_hoogte\hhsk_primaire_kering.shp",
            ),
            (
                r"D:\Work\Project\P1414\GIS\HHSK\Keringen\Regionale_waterkering_1.shp",
                r"D:\Work\Project\P1414\GIS\Keringen_met_hoogte\hhsk_regionale_kering.shp",
            ),
            (
                r"D:\Work\Project\P1414\GIS\HHSK\Keringen\Overige_waterkering_1.shp",
                r"D:\Work\Project\P1414\GIS\Keringen_met_hoogte\hhsk_overige_kering.shp",
            ),
            (
                r"D:\Work\Project\P1414\GIS\WAGV\Keringen\keringen.shp",
                r"D:\Work\Project\P1414\GIS\Keringen_met_hoogte\wagv.shp",
            ),
        ]
    )

    for input_path, output_path in in_out_dict.items():
        gdf = gpd.read_file(input_path).to_crs(crs="EPSG:28992").explode(ignore_index=True)
        print(np.sum(gdf.geometry.type == "MultiLineString"))

        out_gdf = add_height_to_linestrings(gdf=gdf, ahn_path=ahn_path, buffer=11, n_std=1)
        out_gdf.to_file(output_path)
