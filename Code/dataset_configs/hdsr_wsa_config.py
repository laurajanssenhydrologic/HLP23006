from hydrolib.core.io.bc.models import (ForcingBase, ForcingModel,
                                        QuantityUnitPair)
from hydrolib.core.io.ext.models import Boundary, ExtModel, Lateral


class Models:
    class FM:
        one_d_bool = True
        two_d_bool = False
        start_time = 20160601
        stop_time = 86400

        class one_d:
            max_dist_to_struct = 3
            max_snap_dist = 1
            node_distance = 500

        class two_d:
            coupling_type = "2Dto1D"
            dx = 500
            dy = 500
            elevation_raster_path = r"D:\Work\Project\P1414\GIS\AHN\AHN_merged.TIF"
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
                dtmax = 60


class FixedWeirs:
    ## PATHS
    p_folder = r"D:\Work\Project\P1414\GIS"
    flood_defences_path = p_folder + r"\Keringen_met_hoogte\hdsr.shp"

    fixed_weir_index_mapping = dict(
        [
            ("code", "OBJECTID"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
        ]
    )


class RawData:
    ## PATHS
    p_folder = r"D:\Work\Project\HL-23006\GIS"
    # p_folder = r"D:\work\P1414_ROI\GIS"
    branches_path = (
        p_folder + r"\Selectie_watergangen_21062023\sanitized.shp"
    )  # Corrected from V7
    culvert_path = p_folder + r"\Legger\Kokers_Lijnen.shp"
    norm_profile_path = (
        p_folder + r"\Selectie_watergangen_21062023\norm_profiles.shp"
    )

    peil_gebieden_path = p_folder + r"\Legger\BR_Peilgebieden.shp"
    pump_path = p_folder + r"\Legger\Gemalen.shp"
    sluice_path = p_folder + r"\Legger\Sluizen_Lijnen.shp"
    weir_path = dict(
        [
            ("base", p_folder + r"\Legger\BR_Stuwen.shp"),
            # ("concat", p_folder + r"\HDSR\Niet Legger\Kokers_Lijnen_Kerend.shp"),
        ]
    )

    # output_gpkg = p_folder + r"\HDSR\HDSR_hydamo.gpkg"

    ## Branches
    # branch_selection = dict([("column", "CATEGORIEO"), ("value", 1)])
    branch_index_mapping = dict(
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
    )

    ## Bridges
    bridge_index_mapping = dict(
        [
            ("code", "CODE"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("intreeverlies", None),
            ("typeruwheid", None),
            ("ruwheid", None),
            ("uittreeverlies", None),
            ("lengte", "WS_DOORVAA"),
        ]
    )

    ## Culverts
    culvert_index_mapping = dict(
        [
            ("breedteopening", "BREEDTEOPE"),
            ("code", "CODE"),
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
    )

    ## normprofiles
    np_index_mapping = dict(
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
    )
    ## Peil gebied
    peil_index_mapping = dict(
        [
            ("boven peil", ["ZOMERPEIL", "BOVENPEIL"]),
            ("geometry", "geometry"),
            ("onder peil", ["WINTERPEIL", "ONDERPEIL"]),
            ("vast peil", "VASTPEIL"),
        ]
    )

    ## Pumps
    pump_index_mapping = dict(
        [
            ("code", "CODE"),
            ("doelvariabele", None),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("maximalecapaciteit", "MAXIMALECA"),
            ("streefwaarde", None),
            ("peil_marge", None),
        ]
    )

    ## Sluice
    sluice_index_mapping = dict(
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
    )

    ## Weirs
    weir_index_mapping = dict(
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
    )
