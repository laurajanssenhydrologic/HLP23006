from data_structures.config_templates import FMGisFile, RRGisFile
from hydrolib.core.io.bc.models import ForcingBase, ForcingModel, QuantityUnitPair
from hydrolib.core.io.ext.models import Boundary, ExtModel, Lateral

START_TIME = 20010101
STOP_TIME = 86400 * 20
TIMESTEP = 60


class Models:
    class FM:
        one_d_bool = True
        two_d_bool = True
        start_time = START_TIME
        stop_time = STOP_TIME

        class one_d:
            max_dist_to_struct = 3
            max_snap_dist = 1
            node_distance = 500

        class two_d:
            coupling_type = "2Dto1D"
            dx = 100
            dy = 100
            elevation_raster_path = r"D:\Work\Project\P1414\GIS\AHN\AHN_merged.TIF"
            extent_path = r"D:\Work\Project\HL-23006\GIS\Legger\Peilgebieden_dissolved_v2.shp"
            initial_peil_raster_path = r"D:\Work\Project\P1414\GIS\peilen\peilen_jp_25m_full.tif"
            two_d_buffer = 100

        class hydrolib_core_options:
            class external_forcing:
                pass

            class geometry:
                dxmin1d = 1
                usecaching = 1

            class numerics:
                cflmax = 0.7

            class output:
                hisinterval = [0]

            class time:
                dtmax = TIMESTEP

    class RR:
        start_time = START_TIME
        stop_time = STOP_TIME
        timestep = TIMESTEP


P_FOLDER = r"D:\Work\Project\HL-23006\GIS"

duikers = FMGisFile(
    column_mapping=dict(
        [
            ("breedteopening", "BREEDTEOPE"),
            ("code", "CODE"),
            ("doorstroomopening", None),
            ("geometry", "geometry"),
            ("gesloten", None),
            ("globalid", "globalid"),
            ("hoogtebinnenonderkantbene", "HOOGTEBOKB"),
            ("hoogtebinnenonderkantbov", "HOOGTEBO_1"),
            ("hoogteopening", "HOOGTEOPEN"),
            ("intreeverlies", None),
            ("lengte", "LENGTE"),
            ("typeruwheid", None),
            ("ruwheid", None),
            ("uittreeverlies", None),
            ("vormkoker", "VORMKOKER"),
        ]
    ),
    column_selection=dict([("STATUSOBJE", 300)]),
    name="duiker",
    path=P_FOLDER + r"\Legger\Kokers_Lijnen.shp",
)

gemalen = FMGisFile(
    column_mapping=dict(
        [
            ("code", "CODE"),
            ("doelvariabele", None),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("maximalecapaciteit", "MAXIMALECA"),
            ("streefwaarde", None),
            ("peil_marge", None),
        ]
    ),
    column_selection=dict([("STATUSOBJE", 300)]),
    name="gemaal",
    path=P_FOLDER + r"\Legger\Gemalen.shp",
)

keringen = FMGisFile(
    column_mapping=dict([("code", "OBJECTID"), ("geometry", "geometry"), ("globalid", None)]),
    name="keringen",
    path=P_FOLDER + r"\Keringen_met_hoogte\hdsr_simplified.shp",
)

normprofielen = FMGisFile(
    column_mapping=dict(
        [
            ("bodembreedte", "IWS_W_BODB"),
            ("bodemhoogte benedenstrooms", "IWS_W_BODH"),
            ("bodemhoogte bovenstrooms", "IWS_W_BODH"),
            ("code", "CODE"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("hoogte insteek linkerzijde", "IWS_W_INST"),
            ("hoogte insteek rechterzijde", "IWS_W_IN_1"),
            ("taludhelling linkerzijde", "IWS_W_TALU"),
            ("taludhelling rechterzijde", "IWS_W_TA_1"),
            ("typeruwheid", None),
            ("ruwheidhoog", None),
            ("ruwheidlaag", None),
            ("water_width_index", "IWS_W_WATB"),
        ]
    ),
    name="normprofiel",
    path=P_FOLDER + r"\Selectie_watergangen_21062023\norm_profiles.shp",
)

peilgebieden = FMGisFile(
    column_mapping=dict(
        [
            ("boven peil", ["ZOMERPEIL", "BOVENPEIL"]),
            ("geometry", "geometry"),
            ("onder peil", ["WINTERPEIL", "ONDERPEIL"]),
            ("vast peil", "VASTPEIL"),
        ]
    ),
    name="peilgebieden",
    path=P_FOLDER + r"\Legger\BR_Peilgebieden.shp",
)

sluizen = FMGisFile(
    column_mapping=dict(
        [
            ("afvoercoefficient_stuw", None),
            ("afvoercoefficient_opening", None),
            ("code", "CODE"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("hoogstedoorstroombreedte", "DOORVAARTB"),
            ("hoogstedoorstroomhoogte", None),
            ("laagstedoorstroombreedte", "DOORVAARTB"),
            ("laagstedoorstroomhoogte", "KERENDEHOO"),
            ("overlaatonderlaat", None),
            ("soortregelbaarheid", None),
            ("soortstuw", "SOORTSLUIS"),
            ("vormopening", None),
        ]
    ),
    column_selection=dict([("STATUSOBJE", 300)]),
    name="sluis",
    path=P_FOLDER + r"\Legger\Sluizen_Lijnen.shp",
)

stuwen = FMGisFile(
    column_mapping=dict(
        [
            ("afvoercoefficient_stuw", None),
            ("afvoercoefficient_opening", None),
            ("code", "CODE"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("hoogstedoorstroombreedte", "DOORSTROOM"),
            ("hoogstedoorstroomhoogte", "HOOGSTEDOO"),
            ("laagstedoorstroombreedte", "DOORSTROOM"),
            ("laagstedoorstroomhoogte", "LAAGSTEDOO"),
            ("overlaatonderlaat", None),
            ("soortregelbaarheid", "SOORTREGEL"),
            ("soortstuw", "SOORTSTUW"),
            ("vormopening", None),
        ]
    ),
    column_selection=dict([("STATUSOBJE", 300)]),
    name="stuw",
    path=P_FOLDER + r"\Legger\BR_Stuwen.shp",
)

waterlopen = FMGisFile(
    column_mapping=dict(
        [
            ("bodembreedte", "IWS_W_BODB"),
            ("bodemhoogte benedenstrooms", "IWS_W_BODH"),
            ("bodemhoogte bovenstrooms", "IWS_W_BODH"),
            ("code", "CODE"),
            ("geometry", "geometry"),
            ("globalid", "GLOBALID"),
            ("hoogte insteek linkerzijde", "IWS_W_INST"),
            ("hoogte insteek rechterzijde", "IWS_W_IN_1"),
            ("taludhelling linkerzijde", "IWS_W_TALU"),
            ("taludhelling rechterzijde", "IWS_W_TA_1"),
            ("tunnel", False),
            ("typeruwheid", None),
            ("ruwheidhoog", None),
            ("ruwheidlaag", None),
            ("water_width_index", "IWS_W_WATB"),
        ]
    ),
    name="waterloop",
    path=P_FOLDER + r"\Selectie_watergangen_21062023\sanitized.shp",
)
