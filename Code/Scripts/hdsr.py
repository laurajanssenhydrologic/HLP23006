import sys
from pathlib import Path

import geopandas as gpd
import xarray as xr

sys.path.append("D:\Work\git\HL23006\Code")
from data_structures.dhydro_data import DHydroData

folder = r"D:\Work\Project\HL-23006"
ae_path = folder + r"\GIS\koppeling ms sobek\Afwateringseenheden_REFERENTIE_2021-11-12.geojson"
gpkg_file = folder + r"\GIS\HYDAMO\HDSR_v5.gpkg"
gpkg_file_2 = folder + r"\GIS\HYDAMO\HDSR_missing_values.gpkg"
knopen_path = folder + r"\GIS\Riool\comb.shp"
lat_path = folder + r"\GIS\koppeling ms sobek\LATERAL.DAT"
net_path = folder + r"\GIS\koppeling ms sobek\network.nc"
output_folder = folder + r"\Models\HDSR\V07b"
forcing_path = output_folder + r"\dflowfm\boundaryconditions.bc"

config = r"hdsr_wsa_config_v2"
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
    dhd.fm.hydamo_to_gpkg(output_gpkg=gpkg_file)
    dm = dhd.fm.check_missing_values()
    dm.to_gpkg(output_gpkg=gpkg_file_2)
if build_database:
    # 1. initialize an instance of DHydamoData
    dhd = DHydroData()

    # # 2. convert raw data to hydamo data
    dhd.fm.hydamo_from_raw_data(defaults=defaults, config=config, validate=False)
    dhd.fm.clip_structures_by_branches(features=["duiker", "gemaal", "regelmiddel", "stuw"])
    dhd.fm.delete_incomplete_features(
        features=["duiker"], columns=[["hoogtebinnenonderkantbene", "hoogtebinnenonderkantbov"]]
    )
    dhd.fm.validate_dm()

    # 3. save data to gpkg
    dhd.fm.hydamo_to_gpkg(output_gpkg=gpkg_file)

if build_model:
    # 1. initialize an instance of DHydamoData
    dhd = DHydroData()

    # # 2. load data
    dhd.fm.hydamo_from_gpkg(gpkg_file)

    # # 3. save as dhydro model
    dhd.fm.to_dhydro(config=config, output_folder=output_folder)

    dhd.rr.set_basemodel_from_config(config=config)
    dhd.rr.drrmodel.external_forcings.precip = r"D:\Work\Project\HL-23006\Neerslag\default.bui"
    dhd.rr.drrmodel.external_forcings.evap = r"D:\Work\Project\HL-23006\Neerslag\DEFAULT.EVP"

    knopen_gdf = gpd.read_file(knopen_path)

    # extforcefile = dhd.rr.connect_external_weirs_to_branches(
    dhd.rr.connect_external_weirs_to_branches(
        knopen_gdf=knopen_gdf,
        wl_gdf=dhd.fm.dm.waterloop,
        knoop_index_column="objectid_1",
        wl_index_column="code",
    )
    # dhd.fm.merge_ext_files(extforcefile=extforcefile)

    # Load afwateringseenheden and laterals
    ae_gdf = gpd.read_file(ae_path)
    laterals = dhd.rr.read_MF_MS_output(lat_path=lat_path)
    # Save network to temp location and load in xarray
    dhd.fm.fmmodel.geometry.netfile.network.to_file(file=Path(net_path))
    network = xr.open_dataset(net_path)

    # Convert MF-MS input to D-HYDRO
    # extforcefile, forcingmodel = dhd.rr.convert_MF_MS_to_lat
    dhd.rr.convert_MF_MS_to_lat(ae_gdf=ae_gdf, laterals=laterals, network=network)
    # dhd.fm.merge_ext_files(extforcefile=extforcefile)

    dhd.write_dimr(output_folder=output_folder)
    # forcingmodel.save(filepath=forcing_path)
