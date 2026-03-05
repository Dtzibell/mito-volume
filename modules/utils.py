from PyQt6.QtWidgets import QApplication, QFileDialog
import sys
from pathlib import Path
from configparser import ConfigParser


def select_file(type: str):
    config = ConfigParser()
    config.read("config.ini")
    default_dir = config["PATHS"]["DefaultDirectory"]
    app = QApplication(sys.argv)
    # output: [0] list of paths as strings, [1] type of files as string (ex. "*.csv")
    files = QFileDialog.getOpenFileName(
        None,
        f"Select {type} data (.csv's)",
        default_dir,
        "CSV files (*.csv)",
        options = QFileDialog.Option.DontUseNativeDialog,
        )
    return files[0]
