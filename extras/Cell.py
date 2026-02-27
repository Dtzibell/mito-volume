import polars as pl
from polars import col as c
from pathlib import Path

mitocsv = Path(r"/home/tauras/Desktop/Point0000_ChannelmCardinal_Ph-3_Seq0000_s1_run_num1_mCardinal_ref_ch_acdc_output_mask_mitoacdc_outputentation.csv")
ph3csv = Path(r"/home/tauras/Desktop/Point0000_ChannelmCardinal_Ph-3_Seq0000_s1_acdc_output.csv")

class Cell:
    def __init__(self, ph3):
        self.ph3 = ph3
        id = self.getID()
        frame = self.getFrameID()
        self.mito = (pl.scan_csv(mitocsv)
                     .filter((c("Cell_ID") == id) &
                             (c("frame_i") == frame))
                     .collect())
        self.parent_ph3 = (pl.scan_csv(ph3csv)
                       .filter((c("Cell_ID") == self.getParentID()) &
                               (c("frame_i") == self.getFrameID()))
                       .collect())
        self.parent_mito = (pl.scan_csv(mitocsv)
                       .filter((c("Cell_ID") == self.getParentID()) &
                               (c("frame_i") == self.getFrameID()))
                       .collect())

    def getMito(self):
        pass

    def getParentID(self):
        return self.ph3["relative_ID"]

    def getID(self):
        return self.ph3["Cell_ID"]

    def getFrameID(self):
        return self.ph3["frame_i"]

    def getMotherGeneration(self):
        return self.parent_ph3["generation_num"].item()

    def getMitoToVolume(self):
        fluorescence = "mCardinal_concentration_dataPrepBkgr_from_vol_fl_3D"
        volume = "cell_vol_fl"
        self.self_ratio = (self.mito[fluorescence] / self.ph3[volume]).item()

    def getParentMitoToVolume(self):
        fluorescence = "mCardinal_concentration_dataPrepBkgr_from_vol_fl_3D"
        volume = "cell_vol_fl"
        self.parent_ratio = (self.parent_mito[fluorescence] / self.parent_ph3[volume]).item()

    def getSelfToParentRatio(self):
        self.self_to_parent = self.self_ratio / self.parent_ratio
