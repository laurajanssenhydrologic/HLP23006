# import xarray as xr
# import geopandas as gpd
# import pandas as pd


# ds = xr.open_dataset(r"C:\Werk\Projecten\P1414_ROI\Cursus\StandAloneServiceZipfile\Model_database\V20_WBD_v1\dflowfm\network.nc")
# network_x = [ds.network_node_x.values[i] for i in range(len(ds.network_node_x.values))]
# network_y = [ds.network_node_y.values[i] for i in range(len(ds.network_node_y.values))]


# df = pd.DataFrame({'x': network_x , 'y': network_y})
# geometry = gpd.points_from_xy(df.x, df.y, crs="EPSG:28992")

# gdf = gpd.GeoDataFrame(data = {'index': range(len(network_x))}, geometry = geometry)
# print(gdf)

# gdf.to_file(r"C:\Temp\network_x_y_ROI.shp")

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import os
import xarray as xr
# def plot_inundation_volume(nc_path):
#     ncfile = nc.Dataset(nc_path)

#     wd = ncfile.variables['Mesh2d_waterdepth'][:]

#     wd_over_time = np.sum(wd, axis = 1)[:36]

#     # fig, ax = plt.subplots(figsize = (12,12))
#     # ax.plot(
#     #     [x/3 for x in range(len(wd_over_time))],
#     #     wd_over_time
#     #     )
#     # savefigpath = os.path.dirname(os.path.dirname(os.path.dirname(nc_path)))
#     # plt.savefig(os.path.join(savefigpath, 'inundation_total.png'))
#     return np.sum(wd_over_time)

# runs = r"C:\Werk\Projecten\P1414_ROI\Cursus\StandAloneServiceZipfile\Model_runs"
# for folder in os.listdir(runs):
#     map_file = f"{folder}/dflowfm/output/DFM_map.nc"
#     path = os.path.join(runs, map_file)
#     if os.path.exists(path):
#         vol = plot_inundation_volume(path)
#         # print(vol)
#         if vol > 24900:
#             print(path)
#             print(vol)


nc_file = r"C:\Werk\Projecten\P1414_ROI\Cursus-08-05\StandAloneService\Model_database\model_cursus\dflowfm\network.nc"

ds = xr.open_dataset(nc_file)
# print(ds.variables.keys())
# print(ds['mesh1d_node_x'][:])
mesh1d = np.array([(ds.mesh1d_node_x.values[i], ds.mesh1d_node_y.values[i]) for i in range(len(ds.mesh1d_node_x.values))])
mesh2d = np.array([(ds.Mesh2d_face_x.values[i], ds.Mesh2d_face_y.values[i]) for i in range(len(ds.Mesh2d_face_x.values))])

links = ds.links.values


links_geo2 = np.stack((mesh1d[links[:, 0].astype(int)], mesh2d[links[:, 1].astype(int)]), axis=1)

links_geo = []
for link in links: # for each link, draw line between meshd1d point and mesh2d point. 
    new_link = [
        mesh1d[int(link[0])],
        mesh2d[int(link[1])]
    ]
    links_geo.append(new_link)

a = 1

