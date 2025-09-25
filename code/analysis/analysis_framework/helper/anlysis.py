import h5py

class Analyzer:
    """
    Anal eyes
    Analyze
    Anal lies
    """

    def __init__(
        self, log_dir: str = "./", data_dir: str | None = "./", data_type: str = ".h5"
    ) -> None:
        self.log_dir = log_dir
        self.data_dir = data_dir
        self.data_type = data_type

        # ------ Definition of parameters ------
        (
                self.syncSlow,
                self.syncFast2,
                self.lcH,
                self.lcV,
                self.kcH,
                self.syncFast,
                self.sdTrig,
                self.kcV,
                ) = (
                        "ch0",
                        "ch1",
                        "ch2",
                        "ch3",
                        "ch4",
                        "ch5",
                        "ch6",
                        "ch7",
                        )

        self.adt = 0.13  # s - minimum atom trapping duration to be considered "good atom" 

        self.psSave = True  # post selection save
        self.dataSave = True

        # --- Colors ----
        self.colour = {
                "blueDark": (0, 0.3, 0.6),
                "blueLight": (0.5, 0.8, 1),
                "orangeDark": (1, 0.7, 0),
                "orangeLight": (1, 0.8, 0.6),
                "greenDark": (0, 0.6, 0.2),
                "greenLight": (0.7, 1, 0.5),
                "redDark": (0.9, 0, 0),
                "greyLight": (0.7, 0.7, 0.7),
                }

        self.data_dic = None


    def update_data_dir(self, new_data_dir: str) -> None:
        self.data_dir = new_data_dir

    def data_loading(self, path: str | None = None, fileName: str, atom: int):
        """
            INPUTS
            path: directory where to look for the data files
            fileName: names of the data file that will be analysed
            atom: atom number that one wants to analyze
            OUTPUTS
            dataDic: dictionary with the valid timestamps of each channel - dataDic = {'label': [channel,[timeStamp]]}
            """
        if self.data_dir is None:
            print("Please define a data directory first")
            exit(1)

        self.dataDic = {
                "ch0": [0, []],
                "ch1": [1, []],
                "ch2": [2, []],
                "ch3": [3, []],
                "ch4": [4, []],
                "ch5": [5, []],
                "ch6": [6, []],
                "ch7": [7, []],
                }

        if path is None:
            path = self.data_dir

        filedata = h5py.File(path + fileName + ".h5", mode="r")

        for ch in self.dataDic.keys():
            # We put the data into the dictionary
            self.dataDic[ch][1] = (
                    filedata["atom_" + str(atom) + "_" + str(self.dataDic[ch][0])]
                    * filedata.attrs["qu_tau_timebase"]
                    )

        return self.dataDic


