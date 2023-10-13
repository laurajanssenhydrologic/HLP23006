import sys

sys.path.append("D:\Work\git\HL23006\Code")
from data_structures.dhydro_data import DHydroData

folder = r"D:\Work\Project\HL-23006"
gpkg_file = r"D:\Work\Project\HL-23029\GIS\hydroobject_v2.gpkg"
output_folder = folder + r"\Models\HHR\V1"

config = r"hhr_config"
defaults = r"defaults"

validate = False
build_database = False
build_model = True


if build_model:

    # 1. initialize an instance of DHydamoData
    dhd = DHydroData()

    # 2. load data
    dhd.fm.hydamo_from_gpkg(gpkg_file)

    # 3. save as dhydro model
    dhd.fm.to_dhydro(config=config, output_folder=output_folder)
    dhd.write_dimr(output_folder=output_folder)
