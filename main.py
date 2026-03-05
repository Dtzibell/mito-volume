from collections import defaultdict
from time import time
from pathlib import Path
import polars as pl
from polars import col as c
from modules.Cell import Cell
from modules.utils import select_file
from configparser import ConfigParser
import logging
import sys


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
        dct["Size_at_end_of_S"].append(cell.getSizeAtEndOfS())
        dct["G1_length"].append(cell.getLengthOfG1())
        dct["S_length"].append(cell.getLengthOfS())
        dct["Cell_cycle_length"].append(cell.getCellCycleLength())
    else:
        if not cell.parent:
            logger.info(f"Founding mother: ID {int(cell.id)}")
        if not cell.valid:
            logger.info(f"Unfinished cycle: ID {int(cell.id)}")


if __name__ == "__main__":
    logger = logging.getLogger(__name__)
    log_file = "excluded_cells.log"
    logging.basicConfig(
        filename="excluded_cells.log",
        filemode="w",
        encoding="utf-8",
        level=logging.DEBUG,
    )
    config = ConfigParser()
    config.read("config.ini")
    # ph3csv = Path("/home/dtzi/Desktop/Position_0/Images/Point0000_ChannelmCardinal_Ph-3_Seq0000_s1_acdc_output.csv")
    # mitocsv = Path("/home/dtzi/Desktop/Position_0/Images/Point0000_ChannelmCardinal_Ph-3_Seq0000_s1_run_num1_mCardinal_ref_ch_acdc_output_mask_mitoacdc_outputentation.csv")
    ph3csv = Path(select_file("PH3"))
    mitocsv = Path(select_file("Mito"))
    start = time()

    ph3 = pl.scan_csv(ph3csv)
    mito = pl.scan_csv(mitocsv)
    cum_df = (
        mito.join(ph3, on=["Cell_ID", "frame_i"], how="left")
        .select(
            [
                "Cell_ID",
                "frame_i",
                "cell_vol_fl_right",
                "cell_cycle_stage",
                "relative_ID",
                "relationship",
                "mCardinal_concentration_dataPrepBkgr_from_vol_fl_3D",
                "is_history_known",
                "generation_num",
                "time_minutes"
            ]
        )
        .set_sorted("frame_i")
        .collect()
    )
    partitions = cum_df.partition_by("Cell_ID", as_dict=True)

    for p in partitions.values():
        if p.height >= 2:
            imaging_rate = p[1, "time_minutes"] - p[0, "time_minutes"]
            break
        else:
            logging.error("could not determine imaging rate")
            sys.exit(1)

    unknown_history = cum_df.filter(c("is_history_known") == 0).unique(c("Cell_ID"))
    founding_mothers = unknown_history.filter(c("relationship") == "mother").get_column(
        "Cell_ID"
    )

    dct = defaultdict(list)
    for ancestor_id in founding_mothers:
        anc = Cell(ancestor_id, False, partitions, imaging_rate)
        pair_ancestor_to_descendant(anc, dct)

    df = pl.from_dict(dct).sort(c("Bud_ID"))
    df.write_excel(
        Path(config["PATHS"]["OutputDirectory"]) / "output.xlsx",
        freeze_panes=(1, 0),
        autofit=True,
        autofilter=True,
        float_precision=5,
        header_format={"bold": True},
    )
    end = time()
    print(f"Took {round(end - start, 2)} seconds")
    print(f"Excluded cells documented in: {log_file}")
