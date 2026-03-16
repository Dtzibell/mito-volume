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


def pair_ancestor_to_descendants(cell, dct):
    for (i, daughter) in enumerate(cell.daughters):
        pair_ancestor_to_descendants(daughter, dct)
        dct["Bud_ID"].append(int(daughter.id))
        dct["Mother_ID"].append(int(cell.id))
        dct["Mother_generation"].append(i)
        dct["Mito_to_volume_bud"].append(daughter.getMitoToVolumeBud())
        dct["Mito_to_volume_mother"].append(cell.getMitoToVolumeGen(i))
        dct["Bud_to_mother"].append(daughter.getMitoToVolumeBud() / cell.getMitoToVolumeGen(i))
        dct["Size_at_end_of_S"].append(cell.getSizeAtEndOfBud())
        dct["G1_length"].append(cell.getLengthOfG1(i))
        dct["S_length"].append(cell.getLengthOfS(i))
        dct["Cell_cycle_length"].append(cell.getCellCycleLength(i))


if __name__ == "__main__":
    config = ConfigParser()
    config.read("config.ini")
    # ph3csv = Path(
    #     "/home/tauras/Desktop/Point0000_ChannelmCardinal_Ph-3_Seq0000_s1_acdc_output.csv"
    # )
    # mitocsv = Path(
    #     "/home/tauras/Desktop/Point0000_ChannelmCardinal_Ph-3_Seq0000_s1_run_num1_mCardinal_ref_ch_acdc_output_mask_mitoacdc_outputentation.csv"
    # )
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
                "time_minutes",
                "will_divide",
            ]
        )
        .filter(c("will_divide") == 1)
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

    unknown_history_cells = cum_df.filter(c("is_history_known") == 0).unique("Cell_ID")[
        "Cell_ID"
    ]
    founding_mothers = cum_df.filter(
        (c("is_history_known") == 1)
        & (c("relative_ID").is_in(unknown_history_cells.implode()))
    )["Cell_ID"].unique()

    dct = defaultdict(list)
    for id in founding_mothers:
        cell = Cell(id, partitions, imaging_rate)
        pair_ancestor_to_descendants(cell, dct)

    df = pl.from_dict(dct).sort(c("Bud_ID"))
    output_dir = Path(config["PATHS"]["OutputDirectory"])
    output_dir.mkdir(parents=True, exist_ok=True)
    df.write_excel(
        output_dir / "output.xlsx",
        freeze_panes=(1, 0),
        autofit=True,
        autofilter=True,
        float_precision=5,
        header_format={"bold": True},
    )
    end = time()
    print(f"Took {round(end - start, 2)} seconds")
