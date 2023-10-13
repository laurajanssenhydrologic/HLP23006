from pathlib import Path

from hydrolib.dhydamo.io.dimrwriter import DIMRWriter

from data_structures.fm_data import FMData
from data_structures.rr_data import RRData
from data_structures.rtc_data import RTCData

DIMR_BAT_PATH = r"C:\Program Files\Deltares\D-HYDRO Suite 2023.01 1D2D\plugins\DeltaShell.Dimr\kernels\x64\dimr\scripts\run_dimr.bat"


class DHydroData:
    """ """

    def __init__(self):
        self.fm = FMData()
        self.rr = RRData()
        self.rtc = RTCData()

        self.rr.fm = self.fm

    def write_dimr(self, output_folder: str, fm=True, rr=True, rtc=True):
        output_path = Path(output_folder)
        output_path.mkdir(exist_ok=True, parents=True)

        if fm:
            fm = self.fm.write_dimr(output_path=output_path)
        else:
            fm = None

        if rr:
            rr = self.rr.write_drr(output_path=output_path)
        else:
            rr = None

        if rtc:
            rtc = None
        else:
            rtc = None

        dimr = DIMRWriter(dimr_path=DIMR_BAT_PATH, output_path=output_path)
        dimr.write_dimrconfig(fm=fm, rr_model=rr, rtc_model=rtc)
        dimr.write_runbat()
