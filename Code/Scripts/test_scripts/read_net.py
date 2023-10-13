import numpy as np
import xarray as xr
from scipy.spatial import KDTree

net_path = r"D:\Work\Project\HL-23006\models\HDSR\V0\export\dflowfm\network.nc"

ds = xr.open_dataset(net_path)
ids = ds["network_node_id"].values
xs = ds["network_node_x"].values
ys = ds["network_node_y"].values

coords = np.stack([xs, ys], axis=0).T
kdtree = KDTree(data=coords)
dist, nearest = kdtree.query([139003, 444942], k=1)
node_id = ids[nearest].strip()
print(node_id)
