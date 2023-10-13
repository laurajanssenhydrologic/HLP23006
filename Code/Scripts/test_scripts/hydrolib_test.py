from pathlib import Path

from hydrolib.core.basemodel import DiskOnlyFileModel
from hydrolib.core.io.dimr.models import DIMR
from hydrolib.core.io.inifield.models import IniFieldModel, InitialField
from hydrolib.core.io.onedfield.models import OneDFieldGlobal, OneDFieldModel

# assuming the DIMR model is stored in your current working directory.
model_path = Path(r"D:\Work\Project\P1414\Models\ARKNZK\V0\dimr_config.xml")
dimr_model = DIMR(filepath=model_path)
fm_component = dimr_model.component[0]  # Index 0 corresponds with the RRComponent.
fm_model = fm_component.model

f_path = r"D:\Work\Project\P1414\Models\ARKNZK\V00\dflowfm"
fm_model.filepath = f_path + r"\test.mdu"

initialfield = InitialField(
    quantity="waterlevel",
    datafile=DiskOnlyFileModel(filepath=Path("InitialWaterLevel.ini")),
    datafiletype="1dField",
)
initialfieldmodel = IniFieldModel(initial=[initialfield], parameter=[])
# print(fm_model.geometry.inifieldfile.initial)

ifield = fm_model.geometry.inifieldfile.initial[0]
onedfile = OneDFieldModel(global_=OneDFieldGlobal(quantity="WaterLevel", unit="m", value=3))
# onedfile.filepath = f_path + "\\" + fm_model.geometry.inifieldfile.initial[0].datafile
onedfile.filepath = Path(f_path) / fm_model.geometry.inifieldfile.initial[0].datafile.filepath

# onedfile.save()
# print(onedfile.filepath)

# fm_model.save(recurse=True)
# new_path = Path(r"D:\Work\Project\P1414\Models\ARKNZK\V00\dimr_config.xml")
# dimr_model.save(filepath=new_path, recurse=True)
