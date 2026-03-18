import polars as pl
from polars import col as c

# the implementation is not efficient, but it should be fine for now.
class Cell:
    def __init__(
        self, id: int, partitions: dict, IMAGING_RATE: int
    ):
        self.id = id
        # partitions are indexed as tuples (id,)
        self.df = partitions.get((id,), pl.DataFrame())
        self.IMAGING_RATE = IMAGING_RATE
        # further partitions the cell's df by generation number, but it is a list
        # The 0th index of the list is the first generation of the cell - bud
        # stage is excluded from here.
        self.generations = self.df.partition_by("generation_num")[1:]

        # builds the tree:
        # each parent cell nows about its daughter cells, but daughter cells
        # do not know about their parents.
        self.daughters = []
        for gen in self.generations:
            daughter_id = self.getDaughterID(gen)
            self.daughters.append(Cell(daughter_id, partitions, IMAGING_RATE))

    def getDaughterID(self, gen):
        return gen[-1, "relative_ID"]

    def getMitoToVolumeBud(self):
        """
        Retrieves the ratio of mitochondrial fluorescent signal and the cell's
        volume at the end of the cell's budding phase
        """
        fluorescence = "mCardinal_concentration_dataPrepBkgr_from_vol_fl_3D"
        # Apparently, mito datasets have a cell_vol_fl as well. It probably
        # refers to the volume of the mitochondria. After right join of ph3
        # to mito, the ph3 cell_vol_fl is renamed to cell_vol_fl_right
        volume = "cell_vol_fl_right"
        bud_end = self.df.filter(c("generation_num") == 0)[-1]
        #[0, x] retrieves the 0th row, xth column
        return bud_end[0, fluorescence] / bud_end[0, volume]

    def getMitoToVolumeGen(self, gen):
        """
        Retrieves the ratio of mitochondrial fluorescent signal and the cell's
        volume at the end of the cell's generation
        :param gen:
        """
        fluorescence = "mCardinal_concentration_dataPrepBkgr_from_vol_fl_3D"
        # Apparently, mito datasets have a cell_vol_fl as well. It probably
        # refers to the volume of the mitochondria. After right join of ph3
        # to mito, the ph3 cell_vol_fl is renamed to cell_vol_fl_right
        volume = "cell_vol_fl_right"
        gen_end = self.generations[gen][-1] # -1 gets the last row of the generation
        return gen_end[0, fluorescence] / gen_end[0, volume]

    def getSizeAtEndOfBud(self):
        """
        Retrieves the size of the cell at the end of its budding stage
        """
        volume = "cell_vol_fl_right"
        bud_end = self.df.filter(c("generation_num") == 0)[-1]
        return bud_end[0, volume]
    
    def getLengthOfG1(self, gen):
        """
        Retrieves the length of G1
        """
        return (self.generations[gen]
                .filter(c("cell_cycle_stage") == "G1")
                .height * self.IMAGING_RATE)

    def getLengthOfS(self, gen):
        """
        Retrieves the length of S
        """
        return (self.generations[gen]
                .filter(c("cell_cycle_stage") == "S")
                .height * self.IMAGING_RATE
                )

    def getCellCycleLength(self, gen):
        """
        Retrieves the length of the cell cycle.
        """
        return self.getLengthOfG1(gen) + self.getLengthOfS(gen)
