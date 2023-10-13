from data_structures.config_templates import GisFile
from hydrolib.core.io.bc.models import ForcingBase, ForcingModel, QuantityUnitPair
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
