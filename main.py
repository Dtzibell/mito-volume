from collections import defaultdict
from time import time
from pathlib import Path
import polars as pl
from polars import col as c
from extras.Cell import Cell
from extras.utils import select_file
from configparser import ConfigParser

def pair_ancestor_to_descendant(cell: Cell, dct: defaultdict, counter: int = 0):
    for daughter in cell.daughters:
        pair_ancestor_to_descendant(daughter, dct, counter)
    if cell.parent and cell.valid:
        dct["Bud_ID"].append(int(cell.id))
        dct["Mother_ID"].append(int(cell.getParentID()))
        dct["Mother_generation"].append(int(cell.getMotherGeneration()))
        dct["Mito_to_volume_bud"].append(cell.getMitoToVolume())
        dct["Mito_to_volume_mother"].append(cell.getParentMitoToVolume())
        dct["Bud_to_mother"].append(cell.getSelfToParentRatio())

if __name__ == "__main__":
    start = time()
    config = ConfigParser()
    config.read("config.ini")
    ph3csv = Path("/home/dtzi/Desktop/Position_0/Images/Point0000_ChannelmCardinal_Ph-3_Seq0000_s1_acdc_output.csv")
    mitocsv = Path("/home/dtzi/Desktop/Position_0/Images/Point0000_ChannelmCardinal_Ph-3_Seq0000_s1_run_num1_mCardinal_ref_ch_acdc_output_mask_mitoacdc_outputentation.csv")
    # ph3csv = Path(select_file("PH3"))
    # mitocsv = Path(select_file("Mito"))

    ph3 = pl.scan_csv(ph3csv)
    ph3_by_cell = ph3.collect().partition_by("Cell_ID", as_dict=True)
    mito_by_cell = pl.read_csv(mitocsv).partition_by("Cell_ID", as_dict=True)

    cells_ids = ph3.unique("Cell_ID").collect()
    unknown_history = (ph3
                      .filter(c("is_history_known") == 0)
                      .unique(c("Cell_ID"))
                      .collect()
                      )
    founding_mothers = (unknown_history
                      .filter(c("relationship") == "mother")
                      .to_series(1) # intended column Cell_ID
                      )
    
    dct = defaultdict(list)
    for ancestor_id in founding_mothers:
        anc = Cell(ancestor_id, False,
                   ph3_by_cell,
                   mito_by_cell)
        pair_ancestor_to_descendant(anc, dct)

    df = pl.from_dict(dct).sort(c("Bud_ID"))
    df.write_excel("output.xlsx")
    end = time()
    print(f"Took {round(end-start, 2)} seconds")
