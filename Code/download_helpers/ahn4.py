import io
import os.path
import zipfile
from pathlib import Path
from zipfile import BadZipFile

import geopandas as gpd
import numpy as np
import rasterio
import requests
from rasterio.merge import merge
from tqdm import tqdm

# set parameters
ahn3_path = r"D:\work\P1414_ROI\GIS\AHN\AHN3"
ahn4_path = r"D:\work\P1414_ROI\GIS\AHN\AHN4"
ahn_merged_path = r"D:\work\P1414_ROI\GIS\AHN\AHN_merged.TIF"
ahn_tiles_path = r"D:\work\P1414_ROI\GIS\AHN\AHN_bladen.shp"
ahn3_url = r"https://download.pdok.nl/rws/ahn3/v1_0/5m_dtm/"
ahn4_url = r"https://ns_hwh.fundaments.nl/hwh-ahn/ahn4/02b_DTM_5m/"
branches_path = r"D:\work\P1414_ROI\GIS\Randstadmodel_oud\rm_Branches_28992_edited.shp"

max_download_attempts = 2

# manually overide AHN4 in favor of AHN3:
override_list = [
    "32EZ1",
    "32EZ2",
    "32FZ1",
    "32DN2",
    "32GN1",
    "32GN2",
    "32HN1",
    "32DZ1",
    "32DZ2",
    "32GZ1",
    "32GZ2",
    "32HZ1",
    "39AN2",
    "39BN1",
    "39BN2",
    "39EN1",
    "39EN2",
    "39FN1",
    "39AZ2",
    "39BZ1",
    "39BZ2",
    "39EZ1",
    "39EZ2",
    "39FZ1",
]


def ahn_downloader(target_url: str, dst_path: str, max_attempts=2):
    """ """
    try:
        # try to download as many times as max_download_attempts
        for jx in range(max_attempts):
            r = requests.get(target_url)
            if r.ok:
                break

        z = zipfile.ZipFile(io.BytesIO(r.content))
        z.extractall(dst_path)

        return True

    # if file does not exis, return false
    except BadZipFile:
        return False


# create path
Path(ahn3_path).mkdir(exist_ok=True)
Path(ahn4_path).mkdir(exist_ok=True)

# determine extend of AHN to download
branches = gpd.read_file(branches_path)
convex_hull = branches.dissolve(by=None).envelope

# Load ahn tiles and clip with conex_hull of branches
ahn_tiles = gpd.read_file(ahn_tiles_path)
clipped_ahn_tiles = ahn_tiles.clip(convex_hull)


# download all clipped ahn tiles
pbar = tqdm(clipped_ahn_tiles.iterrows(), total=clipped_ahn_tiles.shape[0])
path_list = []
for ix, tile in pbar:
    bladnr = str(tile["bladnr"]).upper()
    ahn3_file_path = ahn3_path + r"\M5_{}.TIF".format(bladnr)
    ahn4_file_path = ahn4_path + r"\M5_{}.TIF".format(bladnr)

    if bladnr in override_list:
        if os.path.isfile(ahn3_file_path):
            pbar.set_description("found AHN3 tile {}".format(bladnr))
            path_list.append(ahn3_file_path)

        else:
            pbar.set_description("downloading AHN3 tile {}".format(bladnr))
            ahn3_tile_url = ahn3_url + r"M5_{}.ZIP".format(bladnr)
            result_ahn3 = ahn_downloader(
                target_url=ahn3_tile_url,
                dst_path=ahn3_path,
                max_attempts=max_download_attempts,
            )

            if result_ahn3:
                path_list.append(ahn3_file_path)
            else:
                pbar.write("AHN3: {} tile not found online")

    else:
        # check if file already exists locally
        if os.path.isfile(ahn4_file_path):
            pbar.set_description("found AHN4 tile {}".format(bladnr))
            path_list.append(ahn4_file_path)

        # else download
        else:
            pbar.set_description("downloading AHN4 tile {}".format(bladnr))
            tile_url = ahn4_url + r"M5_{}.zip".format(bladnr)

            result_ahn4 = ahn_downloader(
                target_url=tile_url, dst_path=ahn4_path, max_attempts=max_download_attempts
            )

            if result_ahn4:
                ahn4_file_path
            else:
                pbar.write("AHN4: {} tile not found online, trying AHN3")
                # continue # prevent adding this path to list

                if os.path.isfile(ahn3_file_path):
                    pbar.set_description("found AHN3 tile {}".format(bladnr))
                    path_list.append(ahn3_file_path)

                ahn3_tile_url = ahn3_url + r"M5_{}.ZIP".format(bladnr)
                result_ahn3 = ahn_downloader(
                    target_url=ahn3_tile_url,
                    dst_path=ahn3_path,
                    max_attempts=max_download_attempts,
                )

                if result_ahn3:
                    path_list.append(ahn3_file_path)
                else:
                    pbar.write("AHN3: {} tile not found online".format(bladnr))


ahn_merged, out_transform = merge(path_list)
ahn_merged[ahn_merged > 100] = np.nan

count, height, width = ahn_merged.shape
profile = {
    "count": count,
    "crs": "EPSG:28992",
    "height": height,
    "width": width,
    "transform": out_transform,
    "dtype": "float32",
}

with rasterio.open(ahn_merged_path, "w", **profile) as dst:
    dst.write(ahn_merged)
