#%%
# # Testfile for overlapping lines
import geopandas as gpd
import shapely
import numpy as np
import pandas as pd
from scipy.spatial import KDTree
import shapelyTools as sT
from importlib import reload
reload(sT)

# Load the watergangen per waterschap
HDSR_water = gpd.read_file(r"D:\work\P1414_ROI\GIS\HYDAMO\HDSR_clipped.gpkg", layer = 'waterloop')
AGV_water = gpd.read_file(r"D:\work\P1414_ROI\GIS\HYDAMO\WAGV_clipped.gpkg", layer = 'waterloop')
HHDL_water = gpd.read_file(r"D:\work\P1414_ROI\GIS\HYDAMO\HHD_clipped.gpkg", layer = 'waterloop')
HHSK_water = gpd.read_file(r"D:\work\P1414_ROI\GIS\HYDAMO\HHSK_clipped.gpkg", layer = 'waterloop')
HHRL_water = gpd.read_file(r"D:\work\P1414_ROI\GIS\HYDAMO\HHR_clipped.gpkg", layer = 'waterloop')

#HDSR_water = gpd.read_file(r"D:\work\P1414_ROI\GIS\Uitgesneden watergangen\HDSR_v4_test.shp")
#AGV_water = gpd.read_file(r"D:\work\P1414_ROI\GIS\Uitgesneden watergangen\AGV_v4_test.shp")
#HHDL_water = gpd.read_file(r"D:\work\P1414_ROI\GIS\Uitgesneden watergangen\HHD_v4_test.shp")
#HHSK_water = gpd.read_file(r"D:\work\P1414_ROI\GIS\Uitgesneden watergangen\HHSK_v4_test.shp")
#HHRL_water = gpd.read_file(r"D:\work\P1414_ROI\GIS\Uitgesneden watergangen\HHR_v4_test.shp")

# Create a list of geometries for both datasets
AGV_geom = [x[1].geometry for x in AGV_water.iterrows()]
HDSR_geom = [x[1].geometry for x in HDSR_water.iterrows()]
HHDL_geom = [x[1].geometry for x in HHDL_water.iterrows()]
HHSK_geom = [x[1].geometry for x in HHSK_water.iterrows()]
HHRL_geom = [x[1].geometry for x in HHRL_water.iterrows()]

agv_mls = shapely.geometry.MultiLineString(AGV_geom)
hdsr_mls = shapely.geometry.MultiLineString(HDSR_geom)
hhdl_mls = shapely.geometry.MultiLineString(HHDL_geom)
hhsk_mls = shapely.geometry.MultiLineString(HHSK_geom)
hhrl_mls = shapely.geometry.MultiLineString(HHRL_geom)

# %%
# Connect all the different datasets at the 10 connection points
# Nr 1: HHDL - HHSK in R'dam Noord, Gordelbrug
hhdl_to_hhsk = sT.snap_endpoints(hhdl_mls, hhsk_mls, max_dist = 5)
print('HHDL en HHSK verbonden in R\'dam Noord')
#%%
# Nr 2: HHDL - HHRL in Noordelijke Sluisbrug Leidschendam
hhdl_to_hhrl_hhsk = sT.snap_endpoints(hhdl_to_hhsk, hhrl_mls, max_dist = 5)
print('HHDL en HHRL verbonden bij Noordelijke Sluisbrug, Leidschendam')
#%%
# Nr 3, 4 & 5: AGV - HDSR op de Vecht, bij Breukelen (ARK) en bij Kockengen
HDSR_buffer = hdsr_mls.buffer(20.0, cap_style=2)
HDSR_buffer_gpd = gpd.GeoDataFrame(HDSR_buffer).rename(columns={0:'geometry'}).set_geometry('geometry')
HDSR_un = HDSR_buffer_gpd.unary_union

AGV_diff = AGV_water.difference(HDSR_un)
AGV_diff_ex = AGV_diff.explode(ignore_index=True)
AGV_diff_ex.drop_duplicates(inplace=True)
AGV_diff_ex_drop = AGV_diff_ex[~AGV_diff_ex.is_empty]
# drop_idx = []
# for i, branch in enumerate(AGV_diff_ex):
#     if branch.length == 0.0:
#         drop_idx.append(i)
# AGV_diff_ex_drop = AGV_diff_ex.drop(drop_idx, axis=0)

agv_to_hdsr = sT.snap_endpoints(list(AGV_diff_ex_drop.geometry), hdsr_mls, max_dist = 80)
print('AGV en HDSR verbonden op de Vecht, bij Breukelen en bij Kockengen')
#%%
# Nr 6: RL- HDSR bij Zwammerdam
RL_buffer = hhrl_mls.buffer(10.0, cap_style=2)
RL_buffer_gpd = gpd.GeoDataFrame(RL_buffer).rename(columns={0:'geometry'}).set_geometry('geometry')
RL_un = RL_buffer_gpd.unary_union

HDSR_diff = HDSR_water.difference(RL_un)
HDSR_diff_ex = HDSR_diff.explode(ignore_index=True)
HDSR_diff_ex.drop_duplicates(inplace=True)
HDSR_diff_ex_drop = HDSR_diff_ex[~HDSR_diff_ex.is_empty]
# drop_idx = []
# for i, branch in enumerate(HDSR_diff_ex):
#     if branch.length == 0.0:
#         drop_idx.append(i)
# HDSR_diff_ex_drop = HDSR_diff_ex.drop(drop_idx, axis=0)
hdsr_to_rl = sT.snap_endpoints(list(HDSR_diff_ex_drop.geometry),hhrl_mls,max_dist = 5)
print('HHRL en HDSR verbonden bij Zwammerdam')

#%%
# Nr 7, 8 & 9: RL - AGV bij de Langeraarsche Plas, Haarlemmervaart en Nieuwe Meer
RL_buffer_15 = hhrl_mls.buffer(15.0, cap_style=2)
RL_buffer_15_gpd = gpd.GeoDataFrame(RL_buffer_15).rename(columns={0:'geometry'}).set_geometry('geometry')
RL_15_un = RL_buffer_15_gpd.unary_union

AGV_to_HDSR_gpd = gpd.GeoDataFrame(agv_to_hdsr).rename(columns={0:'geometry'}).set_geometry('geometry')
AGV_diff_rl = AGV_to_HDSR_gpd.difference(RL_15_un)
AGV_diff_rl_ex = AGV_diff_rl.explode(ignore_index=True)
AGV_diff_rl_ex.drop_duplicates(inplace=True)
AGV_diff_rl_ex_drop = AGV_diff_rl_ex[~AGV_diff_rl_ex.is_empty]
# drop_idx = []
# for i, branch in enumerate(AGV_diff_rl_ex):
#     if branch.length == 0.0:
#         drop_idx.append(i)
# AGV_diff_rl_ex_drop = AGV_diff_rl_ex.drop(drop_idx, axis=0)
agv_to_hdsr_hhrl = sT.snap_endpoints(list(AGV_diff_rl_ex_drop.geometry), hhrl_mls, max_dist=25)
print('HHRL en AGV verbonden bij Langeraarsche Plassen, Haarlemmervaart en Nieuwe Meer')
# %%
all_together = list(hhdl_to_hhrl_hhsk) + list(agv_to_hdsr_hhrl) + \
               list(hdsr_to_rl) + list(hhsk_mls) + list(hhrl_mls)

# %%

all_data_buffered = pd.concat([AGV_water, HDSR_water, HHDL_water, HHRL_water, HHSK_water])
all_data_buffered.geometry = all_data_buffered.geometry.buffer(0.1, cap_style=1)
gdf = gpd.GeoDataFrame(all_together, columns=["geometry"], geometry="geometry", crs=28992)
#%%
# add data from buffered branches if lines in gdf fall within. But keep geometry of branches in gdf
intersected_gdf = gdf.sjoin(all_data_buffered, how="left", predicate="within")
list_of_null = []
for i in intersected_gdf.iterrows():
    if np.isnan(i[1].index_right): # Check whether the original coupling failed
        subgdf = gdf.iloc[[i[0]]]
        
        # Intersected_geo might consists of multiple polygons, especially when connecting at a node
        # Therefore, another selection procedure is needed to get the right branch
        overlay = gpd.overlay(subgdf, all_data_buffered, how='intersection', keep_geom_type=False)
        max_len = 0
        max_idx = 0
        for j in overlay.iterrows():
            if j[1].geometry.length > max_len:
                max_idx = j[0] 
                max_len = j[1].geometry.length
        
        selected_branch_gdf = overlay.iloc[[max_idx]]

        # With then newly selected branches: perform the intersect based sjoin
        intersected_geo = subgdf.sjoin(selected_branch_gdf, how='left',predicate='intersects')
        intersected_gdf = pd.concat([intersected_gdf, intersected_geo])
    
# Remove the NULL entries that caused a problem in the first place
intersected_gdf = intersected_gdf.dropna(axis=0, subset=['index_right'])  
#intersected_gdf_sub = intersected_gdf[['globalid','code', 'geometry', 'typeruwheid']]   
intersected_gdf.to_file(r'D:\work\P1414_ROI\GIS\test_all_intersected.shp')
intersected_gdf.to_file(r'D:\work\P1414_ROI\GIS\test_combined_branches.gpkg', layer="waterloop", driver="GPKG")
print('All data added to the right snapped branches.')
# %%
# Test part to check how this difference function works
HDSR_buffer = hdsr_mls.buffer(20.0, cap_style=2)
HDSR_buffer_gpd = gpd.GeoDataFrame(HDSR_buffer).rename(columns={0:'geometry'}).set_geometry('geometry')
HDSR_un = HDSR_buffer_gpd.unary_union

AGV_diff = AGV_water.difference(HDSR_un)

AGV_diff_ex = AGV_diff.explode(ignore_index=True)
AGV_diff_ex.drop_duplicates(inplace=True)
drop_idx = []
for i, branch in enumerate(AGV_diff_ex):
    if branch.length == 0.0:
        drop_idx.append(i)
AGV_diff_ex_drop = AGV_diff_ex.drop(drop_idx, axis=0)
AGV_diff_ex_drop.to_file('D:\work\P1414_ROI\GIS\AGV_diff_test_v2.shp')
HDSR_buffer_gpd.to_file('D:\work\P1414_ROI\GIS\HDSR_20buf.shp')
HDSR_water.to_file('D:\work\P1414_ROI\GIS\HDSR_test_water.shp')
#agv_to_hdsr = sT.snap_endpoints(AGV_diff, hdsr_mls, max_dist = 80)
# %%
l = [x for x in AGV_diff_ex]
# %%
AGV_test = r"D:\work\P1414_ROI\GIS\Uitgesneden watergangen\AGV_v10_test.shp"
hydrovak = r"D:\work\P1414_ROI\GIS\WAGV\hydrovak.shp"

AGV_gpd = gpd.read_file(AGV_test)
hydrovak_gpd = gpd.read_file(hydrovak)
hydrovak_buffer_gpd = hydrovak_gpd.copy()
hydrovak_buffer_gpd.geometry = hydrovak_buffer_gpd.geometry.buffer(1)
combi = AGV_gpd.sjoin(hydrovak_buffer_gpd, how='left',predicate = 'within')
combi.to_file(r'D:\work\P1414_ROI\GIS\WAGV\Niet legger\hydrovak_combined_v10.shp')
# %%

# Check for duplicate globalids
globlist = []
glob_dup = []
for glob in intersected_gdf.globalid:
    if glob in globlist: 
        glob_dup.append(glob)
    else: 
        globlist.append(glob)
    
print(np.unique(glob_dup, return_counts=True))
# %%
# Check nr of MLS 
westeinder = gpd.read_file(r'D:\work\P1414_ROI\GIS\HHRijnland\Niet legger\Westeinderplassen_cut_v3.shp')
count_mls = 0
mls_idx = []

for n, water in westeinder.iterrows():
    if type(water.geometry) == shapely.geometry.multilinestring.MultiLineString:
        count_mls += 1
        mls_idx.append(n)
# %%
