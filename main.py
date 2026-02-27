from collections import defaultdict
from pathlib import Path
import polars as pl
from polars import col as c
from extras.Cell import Cell
from extras.utils import select_file
from configparser import ConfigParser

if __name__ == "__main__":
    config = ConfigParser()
    config.read("config.ini")
    ph3csv = Path(select_file("PH3"))
    mitocsv = Path(select_file("Mito"))

    ph3 = (pl.scan_csv(ph3csv).reverse() 
           .filter((c("relationship")=="bud") & 
                   (c("frame_i") != 102))
           .unique(c("Cell_ID"))
           .sort(c("Cell_ID"))
           .collect())
    total_cells = ph3.height

    dct = defaultdict(list)

    i = 0
    for row in ph3.iter_rows(named=True):
        if i%20 == 0:
            print(f"Processed {i}/{total_cells} cells")
        cell = Cell(row)
        if cell.mito.height==0:
            print(f"No fluorescence data available for ID {cell.getID()}")
            print(f"Frame of PH3: {cell.getFrameID()}")
            continue
        cell.getMitoToVolume()
        cell.getParentMitoToVolume()
        cell.getSelfToParentRatio()
        dct["Bud_ID"].append(int(cell.getID()))
        dct["Mother_ID"].append(int(cell.getParentID()))
        dct["Mother_generation"].append(int(cell.getMotherGeneration()))
        dct["Mito_to_volume_bud"].append(cell.self_ratio)
        dct["Mito_to_volume_mother"].append(cell.parent_ratio)
        dct["Bud_to_mother"].append(cell.self_to_parent)
        i += 1
    final = pl.from_dict(dct)
    output_dir = Path(config["PATHS"]["OutputDirectory"]) / "output.xlsx"
    final.write_excel(output_dir, 
                      freeze_panes=(1,0),
                      autofit=True,
                      autofilter=True,
                      float_precision=5,
                      header_format={"bold":True}
                      )
