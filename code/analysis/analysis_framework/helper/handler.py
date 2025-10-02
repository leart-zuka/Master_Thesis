import os
from helper.analysis import Analyzer
from typing import List
from rich import print


class AnalysisHandler:
    """
    The main class which will be used in order to anlayze the code at hand
    """

    def __init__(
        self,
        log_dir: str = "./",
        base_data_dir: str = "./",
        year: int = 2025,
        month: int = 1,
        day_topic: str = "",
    ) -> None:
        self.log_dir = log_dir.rstrip("/")
        self.base_data_dir = base_data_dir.rstrip("/")
        self.data_dir: str = self.return_and_load_folder_path(year, month, day_topic)
        self.analzer = Analyzer(log_dir=log_dir, data_dir=self.data_dir)

    def validate_folder_path(self, path_to_folder: str) -> bool:
        path_exists = os.path.exists(path_to_folder)
        if not path_exists:
            print(f"Path [yellow]{path_to_folder}[/yellow] doesn't exist, exiting now")
            exit(1)

        return path_exists

    def return_and_load_folder_path(self, year: int, month: int, day_topic: str) -> str:
        if month < 1 or month > 12:
            print("Please enter a valid month")

        if month < 10:
            month_str = f"0{month}"
        else:
            month_str = str(month)
        path_to_folder = f"{self.base_data_dir}/{year}/{month_str}/{day_topic}"
        self.validate_folder_path(path_to_folder)
        self.data_dir = path_to_folder
        print(f"Loaded up folder [green]{path_to_folder}[/green]")
        return path_to_folder
