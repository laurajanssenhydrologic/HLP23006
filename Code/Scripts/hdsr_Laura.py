#%%
import sys
from pathlib import Path
import os
import shutil
import geopandas as gpd
import xarray as xr

sys.path.append("D:\work\P23006\github\HLP23006\Code")
from data_structures.dhydro_data import DHydroData

folder = r"D:\work\P23006"
gpkg_file_2 = folder + r"\GIS\HYDAMO\HDSR_V8.gpkg"
output_folder = folder + r"\Models\HDSR\V18"

forcing_path = output_folder + r"\dflowfm\boundaryconditions.bc"
boundary_path = folder + r"\GIS\Selectie_watergangen_21062023\boundaries.txt"
waterlevel_path = folder + r"\GIS\peilen\wd_0_v4_test_4_5m.tif"
duiker_path = folder + r"\GIS\Legger\HDSR_Duikersifonhevel.shp"
afsluitmiddel_path = folder + r"\GIS\Legger\HDSR_Afsluitmiddel.shp"
afsluitmiddel_path_out = folder + r"\GIS\Legger\HDSR_Afsluitmiddel_filled.shp"

#rr bestanden
knopen_path = folder + r"\GIS\Riool\comb_v3.shp"
#knopen_path = folder + r"\GIS\Riool\comb_v2_test1.shp"
ae_path = folder + r"\GIS\koppeling ms sobek\Afwateringseenheden_REFERENTIE_2021-11-12.geojson"
lat_path = folder + r"\GIS\koppeling ms sobek\LATERAL.DAT"
net_path = folder + r"\GIS\koppeling ms sobek\network.nc"

config = r"hdsr_wsa_config_vLaura"
defaults = r"defaults"

validate = False
build_database = False
build_model = True

if validate:
    # 1. initialize an instance of DHydamoData
    dhd = DHydroData()

    # # 2. convert raw data to hydamo data
    dhd.fm.hydamo_from_raw_data(defaults=defaults, config=config, validate=False)
    dhd.fm.clip_structures_by_branches()

    # 3. save data to gpkg
    dhd.fm.hydamo_to_gpkg(output_gpkg=gpkg_file_2)
    dm = dhd.fm.check_missing_values()
    dm.to_gpkg(output_gpkg=gpkg_file_2)

if build_database:

    # 1. initialize an instance of DHydamoData
    dhd = DHydroData()

    # # 2. convert raw data to hydamo data
    dhd.fm.fill_afsluitmiddeldata(path_afsluitmiddel=afsluitmiddel_path,path_duiker=duiker_path, path_afsluitmiddel_out=afsluitmiddel_path_out)
    dhd.fm.hydamo_from_raw_data(defaults=defaults, config=config, validate=False)
    dhd.fm.clip_structures_by_branches(features=["duiker", "gemaal", "regelmiddel", "stuw","brug","afsluitmiddel"])
    dhd.fm.delete_incomplete_features(
        features=["duiker"], columns=[["hoogtebinnenonderkantbene", "hoogtebinnenonderkantbov"]])   
    dhd.fm.delete_incomplete_features(
        features=["afsluitmiddel"], columns=[["hoogtebinnenonderkantbene", "hoogtebinnenonderkantbov"]])   
    dhd.fm.validate_dm()

    # 3. save data to gpkg
    dhd.fm.hydamo_to_gpkg(output_gpkg=gpkg_file_2)

if build_model:
    # 1. initialize an instance of DHydamoData
    dhd = DHydroData()

    # # 2. load data
    dhd.fm.hydamo_from_gpkg(gpkg_file_2)

    # # 2.a load boundary conditions for fm part
    forcings_fm, boundaries_fm = dhd.fm.boundary_from_txt(bound_text_path=boundary_path,waterlevel_path=waterlevel_path)

    # # 3. save as dhydro model
    dhd.fm.to_dhydro(config=config, output_folder=output_folder)

    # # 4. add the RR part
    dhd.rr.set_basemodel_from_config(config=config)
    dhd.rr.drrmodel.external_forcings.precip = r"D:\work\P23006\Neerslag\default.bui"
    dhd.rr.drrmodel.external_forcings.evap = r"D:\work\P23006\Neerslag\DEFAULT.EVP"
    knopen_gdf = gpd.read_file(knopen_path)
    dhd.rr.connect_external_weirs_to_branches(knopen_gdf=knopen_gdf,wl_gdf=dhd.fm.dm.waterloop,knoop_index_column="objectid_1",wl_index_column="code")

    # Load afwateringseenheden and laterals
    ae_gdf = gpd.read_file(ae_path)
    laterals = dhd.rr.read_MF_MS_output(lat_path=lat_path)
    # Save network to temp location and load in xarray
    dhd.fm.fmmodel.geometry.netfile.network.to_file(file=Path(net_path))
    network = xr.open_dataset(net_path)

    # Convert MF-MS input to D-HYDRO
    dhd.rr.convert_MF_MS_to_lat(ae_gdf=ae_gdf, 
                                laterals=laterals, 
                                network=network,
                                boundary_fm = boundaries_fm,  
                                forcing_fm = forcings_fm,      
                                forcing_path=output_folder + r"\dflowfm\boundaryconditions_v2.bc")
    #dhd.fm.merge_ext_files(extforcefile=extforcefile)

    dhd.write_dimr(output_folder=output_folder,rr=True)

    # Fix boundary conditions filename
    os.remove(forcing_path)
    os.rename(output_folder + r"\dflowfm\boundaryconditions_v2.bc",forcing_path)
    shutil.copyfile(folder+r"\GIS\default.tmp",output_folder+r"\rr\default.tmp")
    # forcingmodel.save(filepath=forcing_path)

# %%
