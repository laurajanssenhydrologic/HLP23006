import geopandas as gpd
import numpy as np
import rasterio

duikers_path = r"D:\Work\Project\HL-23006\GIS\Legger\Kokers_Lijnen.shp"
duikers_gdf = gpd.read_file(duikers_path).to_crs("EPSG:28992")
ahn_path = r"D:\Work\Project\P1414\GIS\AHN\AHN4_WSS_filled.tif"

ahn = rasterio.open(ahn_path)

duikers_gdf.loc[duikers_gdf["HOOGTEBOKB"] < -5, "HOOGTEBOKB"] = np.nan
duikers_gdf.loc[duikers_gdf["HOOGTEBO_1"] < -5, "HOOGTEBO_1"] = np.nan
duikers_gdf.loc[duikers_gdf["HOOGTEBOKB"] > 25, "HOOGTEBOKB"] = np.nan
duikers_gdf.loc[duikers_gdf["HOOGTEBO_1"] > 25, "HOOGTEBO_1"] = np.nan

duikers_gdf["HOOGTEBOK1"] = duikers_gdf["HOOGTEBOKB"]
duikers_gdf["HOOGTEBOK2"] = duikers_gdf["HOOGTEBO_1"]
for ix, duiker in duikers_gdf.iterrows():
    nanmean = np.nanmean([duiker["HOOGTEBOKB"], duiker["HOOGTEBO_1"]])
    if nanmean == 0:
        nanmean = np.nan
    # if nanmean != 0:
    #     duikers_gdf.at[ix, "HOOGTEBOK"] = nanmean

    if np.isnan(duiker["HOOGTEBOK1"]):
        duikers_gdf.at[ix, "HOOGTEBOK1"] = nanmean
    if np.isnan(duiker["HOOGTEBOK2"]):
        duikers_gdf.at[ix, "HOOGTEBOK2"] = nanmean
    x, y = duiker.geometry.xy

    duikers_gdf.at[ix, "X1"] = x[-1]
    duikers_gdf.at[ix, "Y1"] = y[-1]
    duikers_gdf.at[ix, "AHN1"] = next(ahn.sample([(x[-1], y[-1])]))
    duikers_gdf.at[ix, "X2"] = x[0]
    duikers_gdf.at[ix, "Y2"] = y[0]
    duikers_gdf.at[ix, "AHN2"] = next(ahn.sample([(x[0], y[0])]))

ahn.close()
duikers_gdf.to_file(r"D:\Work\Project\HL-23006\GIS\Legger\Kokers_Lijnen_filled.shp")
