from __future__ import annotations
import polars as pl
from polars import col as c

class Cell:

    def __init__(self, id: int, parent: Cell | bool,
                 ph3: pl.LazyFrame, mito: pl.LazyFrame):
        self.id = id
        self.ph3 = (ph3
                    .filter(c("Cell_ID") == id)
                    .collect()
                    )
        self.mito = (mito
                    .filter(c("Cell_ID") == id)
                    .collect()
                     )
        self.budding = self.ph3.filter(c("cell_cycle_stage") == "S")

        # sets valid flag for cells that have finished a full S phase.
        if self.budding.height < self.ph3.height:
            self.valid = True
        else: self.valid = False
        if self.valid and parent:
            self.ph3_bud_end = (self.budding
                                     .filter(c("relationship") == "bud")
                                     .reverse()
                                     .unique(c("Cell_ID"))
                                     )
            self.mito_bud_end = (self.mito
                                 .filter(c("frame_i") == 
                                         self.ph3_bud_end[0, "frame_i"])
                                 )

        daughters = (self.ph3
                     .filter((c("relationship") == "mother") & (c("cell_cycle_stage") == "S"))
                     .unique(c("relative_ID"))
                     .to_series(6) # intended column Cell_ID
                     )
        self.parent = parent
        self.daughters = []
        for daughter_id in daughters:
            self.daughters.append(Cell(daughter_id, self, ph3, mito))

    def getParentID(self):
        match self.parent:
            case Cell() as parent:
                return parent.id
            case bool():
                print(f"Cell {self.id} has no ancestor")
                return -1

    def getBudEndFrame(self):
        return self.ph3_bud_end[0, "frame_i"]

    def getMotherGeneration(self):
        match self.parent:
            case Cell() as parent:
                return (parent.ph3
                        .filter(c("frame_i") == self.getBudEndFrame())
                        )[0, "generation_num"]
            case bool():
                print(f"Cell {self.id} has no ancestor")
                return -1

    def getMitoToVolume(self):
        fluorescence = "mCardinal_concentration_dataPrepBkgr_from_vol_fl_3D"
        volume = "cell_vol_fl"
        self.self_ratio = (self.mito_bud_end[0,fluorescence] /
                           self.ph3_bud_end[0,volume])
        return self.self_ratio

    def getParentMitoToVolume(self):
        fluorescence = "mCardinal_concentration_dataPrepBkgr_from_vol_fl_3D"
        volume = "cell_vol_fl"
        match self.parent:
            case Cell() as parent:
                parent_mito = (parent.mito
                               .filter(c("frame_i") == self.getBudEndFrame())
                               )
                parent_ph3 = (parent.ph3
                              .filter(c("frame_i") == self.getBudEndFrame())
                              )
            case bool():
                print(f"Cell {self.id} has no ancestor")
                self.parent_ratio = -1
                return
        self.parent_ratio = parent_mito[0, fluorescence] / parent_ph3[0, volume]
        return self.parent_ratio

    def getSelfToParentRatio(self):
        self.self_to_parent = self.self_ratio / self.parent_ratio
        return self.self_to_parent


