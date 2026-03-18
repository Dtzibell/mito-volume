from collections import defaultdict
from time import time
from pathlib import Path
import polars as pl
from polars import col as c
from modules.Cell import Cell
from modules.utils import select_file
from configparser import ConfigParser
import sys

def pair_ancestor_to_descendants(cell, dct):
    """
    Retrieves data about each cell cycle that relates to the given cell.
    Recursively goes down the cell's genealogical tree.

    :param cell: a Cell object
    :param dct: a dictionary where the cell cycles are to be saved
    """
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

def parse_config(config_file_name: str) -> ConfigParser:
    """
    Retrieves the config

    :param config_file_name: name of the file where the configuration is stored.
    """
    config = ConfigParser()
    config.read(config_file_name)
    print(config.sections())
    return config

if __name__ == "__main__":
    config = parse_config("config.ini")
    # convenience for development
    # ph3csv = Path("/home/tauras/Desktop/Point0000_ChannelmCardinal_Ph-3_Seq0000_s1_acdc_output.csv")
    # mitocsv = Path("/home/tauras/Desktop/Point0000_ChannelmCardinal_Ph-3_Seq0000_s1_run_num1_mCardinal_ref_ch_acdc_output_mask_mitoacdc_outputentation.csv")
    
    # file selection windows. Can modify to automatically select _output and 
    # _outputentation via pathlib.Path.walk(), but need to be sure that names
    # are always that. Hence not implemented.
    ph3csv = Path(select_file("PH3"))
    mitocsv = Path(select_file("Mito"))

    # starts run timer, I like it for performance evaluation
    start = time()

    # scan_csv creates polars.LazyFrames, which are not real DataFrames.
    # They only become real DataFrames once collect is called
    ph3 = pl.scan_csv(ph3csv)
    mito = pl.scan_csv(mitocsv)

    # ph3 dataset is merged on the mito dataset. "left" indicates that mito is 
    # the reference dataset - only rows that exist in mito will appear in the 
    # cumulative dataframe (hence the problem with missing first frame for buds)
    cum_lf = (mito.join(ph3, on=["Cell_ID", "frame_i"], how="left")
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
              )
    # It is possible to merge on ph3, but it causes many problems and is super
    # dirty... Do not recommend.

    # removes the cells that do not divide
    cum_df = (cum_lf.filter(c("will_divide") == 1)
              .collect() # makes the DataFrame
              )

    # splits the dataframe into a dictionary of dataframes for each cell_id
    partitions = cum_df.partition_by("Cell_ID", as_dict=True)

    # finds the experiments imaging rate
    for p in partitions.values():
        if p.height >= 2:
            imaging_rate = p[1, "time_minutes"] - p[0, "time_minutes"]
            break
        else:
            print("could not determine imaging rate")
            sys.exit(1)

    # founding mothers are cells that do not have progenitors with a known
    # history. For example, if cell id 1 does not have a known history, all
    # of its buds throughout the experiment are labeled founding mothers.
    unknown_history_cells = cum_df.filter(c("is_history_known") == 0).unique("Cell_ID")[
        "Cell_ID"
    ]
    founding_mothers = cum_df.filter(
        (c("is_history_known") == 1)
        # dont know why implode has to be called here, but theres this
        # https://github.com/pola-rs/polars/pull/22178
        & (c("relative_ID").is_in(unknown_history_cells.implode()))
    )["Cell_ID"].unique()

    # determining founding_mothers allows me to build a tree. E.g.: 
    #            1 <-- founding mother
    #           / \
    #          3  2
    #         /  / \ ... other cells
    # working with trees is very intuitive, hence the structure.

    # defaultdict(list) allows list operations to be called on each item regardless
    # of whether the item is present or not.
    dct = defaultdict(list)

    for id in founding_mothers:
        # refer to modules/Cell.py for documentation
        cell = Cell(id, partitions, imaging_rate)
        pair_ancestor_to_descendants(cell, dct)

    # creates the output csv and terminates the script
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
