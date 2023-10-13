import os
import sys
from typing import Optional

sys.path.append("D:\Work\git\GIS_tools\Code")

import pandera as pa
from data_structures.dhydamo_data_model_checks import (
    geometry_check,
    globalid_check,
    none_geometry_check,
    validate_codes,
)
from data_structures.hydamo_globals import (
    HYDAMO_SHAPE_NUMS,
    MANAGEMENT_DEVICE_TYPES,
    ROUGHNESS_MAPPING_LIST,
)
from pandera.typing import DataFrame, Series
from pandera.typing.geopandas import GeoSeries
from pydantic import BaseModel, root_validator


class BasicSchema(pa.SchemaModel):
    class Config:
        coerce = True  # allow datatype conversion
        name = "BasicSchema"
        strict = True  # no additional columns allowed


class PDBasicShema(BasicSchema):
    globalid: Series[str] = pa.Field(unique=True)
    geometry: GeoSeries = pa.Field(nullable=True)

    @pa.check("globalid", element_wise=True, name="globalidcheck")
    def _globalid_check(cls, globalid: str) -> bool:
        return globalid_check(globalid=globalid)

    @pa.check("geometry", element_wise=True, name="None geometry check")
    def _none_geometry_check(cls, geometry: None) -> bool:
        return none_geometry_check(geometry=geometry)


class GPDBasicShema(BasicSchema):
    globalid: Series[str] = pa.Field(unique=True)
    geometry: GeoSeries

    @pa.check("globalid", element_wise=True, name="globalidcheck")
    def _globalid_check(cls, globalid: str) -> bool:
        return globalid_check(globalid=globalid)

    @pa.check("geometry", element_wise=True, name="Point or line geometry check")
    def _geometry_check(cls, geometry: None) -> bool:
        return geometry_check(geometry=geometry)


class BrugSchema(GPDBasicShema):
    code: Optional[Series[str]] = pa.Field(unique=True)
    intreeverlies: Series[float] = pa.Field(ge=0, le=1)
    lengte: Series[float] = pa.Field(gt=0)
    ruwheid: Series[float] = pa.Field(gt=0)
    typeruwheid: Series[str] = pa.Field(isin=ROUGHNESS_MAPPING_LIST)
    uittreeverlies: Series[float] = pa.Field(ge=0, le=1)


class HydamoDataModel(BaseModel):
    brug: Optional[DataFrame[BrugSchema]]


bs = BrugSchema
print(bs.to_schema().columns)
print(list(bs.to_schema().columns))
print(bs.to_schema())

# print(bs.to_yaml())
