import os
from helper.anlysis import Analyzer


class Analysis_Handler:
    """
    The main class which will be used in order to anlayze the code at hand
    """

    def __init__(self, log_dir: str = "./", base_data_dir: str = "./") -> None:
        self.log_dir = log_dir
        self.base_data_dir = base_data_dir
        self.data_dir: str | None = None
        self.analzer = Analyzer(log_dir=log_dir)

    def validate_folder_path(self, path_to_folder: str) -> bool:
        path_exists = os.path.exists(path_to_folder)
        if not path_exists:
            print("Path doesn't exist, exiting now")
            exit(1)
        return path_exists

    def return_and_load_folder_path(self, year: int, month: int, day_topic: str) -> str:
        if month < 10:
            month_str = f"0{month}"
        else:
            month_str = str(month)
        path_to_folder = f"{self.base_data_dir}/{year}/{month_str}/{day_topic}/"
        self.validate_folder_path(path_to_folder)
        self.data_dir = path_to_folder
        self.analzer.update_data_dir(self.data_dir)
        return path_to_folder
