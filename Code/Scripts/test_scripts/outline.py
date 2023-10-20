import geopandas as gpd
import matplotlib.pyplot as plt

in_path = r"D:\Work\Project\HL-23006\GIS\Legger\BR_Peilgebieden.shp"
out_path = r"D:\Work\Project\HL-23006\GIS\Legger\Peilgebieden_dissolved_v2.shp"

gdf = gpd.read_file(in_path)
gdf = gdf.to_crs("EPSG:28992")

gdf["geometry"] = gdf["geometry"].buffer(0.1)
gdf = gdf.dissolve(by=None)
gdf.to_file(out_path)
gdf.plot()
plt.show()
