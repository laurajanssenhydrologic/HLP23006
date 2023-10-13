# %% Init
print("Initializing")

import geopandas as gpd

BUFFER_DIST = 50
EPSG = 28992


def clip_profiles(in_profiles_path: str, buffered_branches: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    in_branches = gpd.read_file(in_profiles_path).to_crs(buffered_branches.crs)
    return gpd.overlay(in_branches, buffered_branches, how="intersection", keep_geom_type=True)


def read_rm_branches(rm_branches_path: str) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    rm_branches = gpd.read_file(rm_branches_path).to_crs(crs=EPSG)
    ondd_bool = (rm_branches["Source"].str.contains("OND")) | (
        rm_branches["Target"].str.contains("OND")
    )
    onderdoorgangen = rm_branches.loc[ondd_bool, :]
    rm_branches = rm_branches.loc[~ondd_bool, :]
    return rm_branches, onderdoorgangen


old_rm_branches_path = r"D:\work\P1414_ROI\GIS\Randstadmodel_oud\rm_Branches_28992.shp"

old_rm_branches, onderdoorgangen = read_rm_branches(old_rm_branches_path)
buffered_old_rm_branches = gpd.GeoDataFrame(geometry=old_rm_branches.buffer(distance=BUFFER_DIST))


# %% AGV
print("AGV")

agv_profiles_path = (
    r"D:\work\P1414_ROI\GIS\WAGV\metingprofielpunt_v13\metingprofielpunt_v13_clipped.shp"
)
clipped_agv_profiles_path = (
    r"D:\work\P1414_ROI\GIS\WAGV\metingprofielpunt_v13\metingprofielpunt_v13_clipped_rm.shp"
)

intersected_profiles = clip_profiles(
    in_profiles_path=agv_profiles_path, buffered_branches=buffered_old_rm_branches
)
intersected_profiles.to_file(clipped_agv_profiles_path)
