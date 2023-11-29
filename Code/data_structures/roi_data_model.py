import os
from typing import Any, Optional

import pandera as pa
from pandera.typing import DataFrame, Series
from pydantic import BaseModel, root_validator

from data_structures import hydamo_data_model as hdm
from data_structures.dhydamo_data_model_checks import validate_codes
from data_structures.hydamo_globals import HYDAMO_WEIR_TYPES, ROUGHNESS_MAPPING_LIST


class BrugSchema(hdm.BrugSchema):
    pass


class DuikerSchema(hdm.DuikerSchema):
    pass


class GemaalSchema(hdm.GemaalSchema):
    pass


class Hydroobject_normgpSchema(hdm.Hydroobject_normgpSchema):
    pass


class KunstwerkopeningSchema(hdm.KunstwerkopeningSchema):
    pass


class NormgeparamprofielwaardeSchema(hdm.NormgeparamprofielwaardeSchema):
    pass


class PompSchema(hdm.PompSchema):
    code: Series[str] = pa.Field(unique=True)  # addition to confluence


class ProfielgroepSchema(hdm.ProfielgroepSchema):
    pass


class ProfiellijnSchema(hdm.ProfiellijnSchema):
    pass


class ProfielpuntSchema(hdm.ProfielpuntSchema):
    pass


class RegelmiddelSchema(hdm.RegelmiddelSchema):
    code: Series[str]  # addition to confluence


class RuwheidsprofielSchema(hdm.RuwheidsprofielSchema):
    code: Series[str]  # addition to confluence


class SturingSchema(hdm.SturingSchema):
    code: Series[str] = pa.Field(unique=True)  # addition to confluence
    doelvariabele: Series[str]  # addition to confluence


class StuwSchema(hdm.StuwSchema):
    soortstuw: Series[float] = pa.Field(isin=HYDAMO_WEIR_TYPES)  # addition to confluence


class WaterloopSchema(hdm.WaterloopSchema):
    peil: Optional[Series[float]] = pa.Field(nullable=True)
    tunnel: Optional[Series[bool]] = pa.Field(nullable=True)
    typeruwheid: Series[str] = pa.Field(isin=ROUGHNESS_MAPPING_LIST)  # addition to confluence

    class Config:
        strict = "filter"


class FMDataModel(BaseModel):
    """
    Datamodel that validates attributes using pydantic and pandera, using pandera SchemaModels.
    The former only checks for type (i.e. DataFrame), the latter checks the content of the DataFrames
    Ensures that data is compatible with DHydamo

    Attributes:
        brug (gpd.GeoDataFrame): geodataframe containing bridge data
        duiker (gpd.GeoDataFrame): geodataframe containing culvert data
        gemaal (gpd.GeoDataFrame): geodataframe containing pump station data
        hydroobject_normgp (gpd.GeoDataFrame): geodataframe containing norm profile id data
        kunstwerkopening (gpd.GeoDataFrame): geodataframe containing weir opening data
        normgeparamprofielwaarde (gpd.GeoDataFrame): geodataframe containing norm profile data
        pomp (gpd.GeoDataFrame): geodataframe containing pump data
        profielgroep (gpd.GeoDataFrame): geodataframe containing measured profile groups
        profiellijn (gpd.GeoDataFrame): geodataframe containing measured profile lines
        profielpunt (gpd.GeoDataFrame): geodataframe containing measured profile points
        regelmiddel (gpd.GeoDataFrame): geodataframe containing management device data for weirs
        ruwheidsprofiel (gpd.GeoDataFrame): geodataframe containing roughness data for measured profiles
        sturing (gpd.GeoDataFrame): geodataframe containing management data for pumps
        stuw (gpd.GeoDataFrame): geodataframe containing weir data
        waterloop (gpd.GeoDataFrame): geodataframe containing branch data
    """

    brug: Optional[DataFrame[BrugSchema]]
    duiker: Optional[DataFrame[DuikerSchema]]
    afsluitmiddel: Optional[DataFrame[DuikerSchema]]
    doorbraak: Optional[Any]
    gemaal: Optional[DataFrame[GemaalSchema]]
    hydroobject_normgp: Optional[DataFrame[Hydroobject_normgpSchema]]
    keringen: Optional[Any]
    kunstwerkopening: Optional[DataFrame[KunstwerkopeningSchema]]
    normgeparamprofielwaarde: Optional[DataFrame[NormgeparamprofielwaardeSchema]]
    overiglijnelement: Optional[Any]
    peilgebieden: Optional[Any]
    pomp: Optional[DataFrame[PompSchema]]
    profielgroep: Optional[DataFrame[ProfielgroepSchema]]
    profiellijn: Optional[DataFrame[ProfiellijnSchema]]
    profielpunt: Optional[DataFrame[ProfielpuntSchema]]
    regelmiddel: Optional[DataFrame[RegelmiddelSchema]]
    rivier_profielen: Optional[Any]
    rivier_profielen_data: Optional[Any]
    rivier_ruwheid: Optional[Any]
    ruwheidsprofiel: Optional[DataFrame[RuwheidsprofielSchema]]
    sturing: Optional[DataFrame[SturingSchema]]
    stuw: Optional[DataFrame[StuwSchema]]
    waterloop: Optional[DataFrame[WaterloopSchema]]

    @root_validator(pre=False)
    def validate_codes(cls, values):
        """
        Validate that codes are unique amongst structures
        """

        return validate_codes(values)

    class Config:
        extra = "forbid"  # do not allow non-specified fields to be assigned
        validate_assignment = True  # validates new fields upon assignment

    def to_gpkg(self, output_gpkg: str) -> None:
        """
        Class method that saves fields that are not None to gpkg

        Args:
            output_gpkg (str): file location to write geopackge to

        Returns:
            None
        """
        try:
            os.remove(output_gpkg)
        except FileNotFoundError:
            pass
        except OSError:
            print("could not acces gpkg")

        for key, value in self.__dict__.items():
            if value is not None:
                gdf = getattr(self, key)
                if gdf.index.name in gdf.columns:
                    gdf.reset_index(drop=True, inplace=True)

                try:
                    gdf.to_file(output_gpkg, layer=key, driver="GPKG")
                except Exception as e:
                    print("\nfailed to save {}\n".format(key))
                    print(e)
