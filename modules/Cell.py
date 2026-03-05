from __future__ import annotations
import polars as pl
from polars import col as c
from polars.exceptions import ColumnNotFoundError


class Cell:
    def __init__(
        self, id: int, parent: Cell | bool, partitions: dict, IMAGING_RATE: int
    ):
        self.id = id
        self.df = partitions.get((id,), pl.DataFrame())
        self.IMAGING_RATE = IMAGING_RATE
        self.parent = parent

        self.first_cell_cycle = self.df.filter(c("generation_num").is_in([0, 1]))
        self.budding = self.first_cell_cycle.filter(c("cell_cycle_stage") == "S")

        # sets valid flag for cells that have finished a full S phase.
        if self.budding.height < self.df.height:
            self.valid = True
        else:
            self.valid = False
        if self.valid and parent:
            self.bud_end = self.budding.filter(c("relationship") == "bud")[-1]

        daughters = (
            self.df.lazy()
            .filter((c("relationship") == "mother") & (c("cell_cycle_stage") == "S"))
            .unique(c("relative_ID"))
            .collect()
            .get_column("relative_ID")
        )
        self.daughters = []
        for daughter_id in daughters:
            try:
                self.daughters.append(Cell(daughter_id, self, partitions, IMAGING_RATE))
            except ColumnNotFoundError:
                print(
                    f"{int(daughter_id)} is not present within the mitochondria dataset"
                )
                continue

    def getParentID(self):
        match self.parent:
            case Cell() as parent:
                return parent.id
            case bool():
                print(f"Cell {self.id} has no ancestor")
                return -1

    def getBudEndFrame(self):
        return self.bud_end[0, "frame_i"]

    def getMotherGeneration(self):
        match self.parent:
            case Cell() as parent:
                return (parent.df.filter(c("frame_i") == self.getBudEndFrame()))[
                    0, "generation_num"
                ]
            case bool():
                print(f"Cell {self.id} has no ancestor")
                return -1

    def getMitoToVolume(self):
        fluorescence = "mCardinal_concentration_dataPrepBkgr_from_vol_fl_3D"
        # Apparently, mito datasets have a cell_vol_fl as well. It probably
        # refers to the volume of the mitochondria. After right join of ph3
        # to mito, the ph3 cell_vol_fl is renamed to cell_vol_fl_right
        volume = "cell_vol_fl_right"

        self.self_ratio = self.bud_end[0, fluorescence] / self.bud_end[0, volume]
        return self.self_ratio

    def getParentMitoToVolume(self):
        fluorescence = "mCardinal_concentration_dataPrepBkgr_from_vol_fl_3D"
        # Apparently, mito datasets have a cell_vol_fl as well. It probably
        # refers to the volume of the mitochondria. After right join of ph3
        # to mito, the ph3 cell_vol_fl is renamed to cell_vol_fl_right
        volume = "cell_vol_fl_right"
        match self.parent:
            case Cell() as parent:
                parent_frame = parent.df.filter(c("frame_i") == self.getBudEndFrame())
            case bool():
                print(f"Cell {self.id} has no ancestor")
                self.parent_ratio = -1
                return
        self.parent_ratio = parent_frame[0, fluorescence] / parent_frame[0, volume]
        return self.parent_ratio

    def getSelfToParentRatio(self):
        self.self_to_parent = self.self_ratio / self.parent_ratio
        return self.self_to_parent

    def getLengthOfS(self):
        return self.budding.height * self.IMAGING_RATE

    def getLengthOfG1(self):
        return (self.first_cell_cycle.height - self.budding.height) * self.IMAGING_RATE

    def getCellCycleLength(self):
        return self.first_cell_cycle.height * self.IMAGING_RATE

    def getSizeAtEndOfS(self):
        return self.budding[-1, "cell_vol_fl_right"]
