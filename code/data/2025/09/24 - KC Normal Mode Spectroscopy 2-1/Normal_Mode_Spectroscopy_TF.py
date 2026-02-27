import os
import sys
import time
import h5py
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np
import lmfit
from scipy.optimize import curve_fit, minimize
import pandas as pd
import pickle
import matplotlib.colors as mcolors
from typing import List


class DataLoading:
    def __init__(self):
        # ------ Definition of parameters ------
        print("Data Loading")

    # ------ Data loading ------
    def data_loading(self, path: str, fileName: str, atom: int):
        """
        INPUTS
        path: directory where to look for the data files
        fileName: names of the data file that will be analysed
        atom: atom number that one wants to analyze
        OUTPUTS
        dataDic: dictionary with the valid timestamps of each channel - dataDic = {'label': [channel,[timeStamp]]}
        """
        dataDic = {
            "ch0": [0, []],
            "ch1": [1, []],
            "ch2": [2, []],
            "ch3": [3, []],
            "ch4": [4, []],
            "ch5": [5, []],
            "ch6": [6, []],
            "ch7": [7, []],
        }

        filedata = h5py.File(path + "/" + fileName + ".h5", mode="r")

        for ch in dataDic.keys():
            # We put the data into the dictionary
            dataDic[ch][1] = (
                filedata["atom_" + str(atom) + "_" + str(dataDic[ch][0])]
                * filedata.attrs["qu_tau_timebase"]
            )

        filedata.close()

        return dataDic

    # ------ Data loading from the good atom periods ------
    def data_goodAtoms(self, path: str, fileName: str, atomDic: dict):
        """
        INPUTS
        path: directory where to look for the data files
        fileName: names of the data file that will be analysed
        atomDic: dictionary with the "good" atoms and the corresponding time period
        OUTPUTS
        dataVD: dictionary with the valid timestamps of each channel - dataDic = {'label': [channel,[timeStamp]]}
        """
        dataVD = {
            "ch0": [0, []],
            "ch1": [1, []],
            "ch2": [2, []],
            "ch3": [3, []],
            "ch4": [4, []],
            "ch5": [5, []],
            "ch6": [6, []],
            "ch7": [7, []],
        }

        for atom in tqdm(atomDic.keys(), file=sys.stdout):
            (dataD) = self.data_loading(path, fileName, atom)

            dataD["ch1"][1] = np.array([x for x in dataD["ch1"][1]])

            # We take the time-stamps from the valid atom time period
            right_sf = np.searchsorted(
                dataD["ch5"][1], dataD["ch0"][1][0] + atomDic[atom][1]
            )  # index of the first syncFastScan that we don't consider
            left_sf = np.searchsorted(
                dataD["ch5"][1], dataD["ch0"][1][0] + atomDic[atom][0]
            )  # index of the first syncFastScan that we consider

            timeInit = (
                dataD["ch5"][1][left_sf] - 1e-9
            )  # time after which we start to consider data
            timeEnd = (
                dataD["ch5"][1][right_sf] - 1e-9
            )  # time after which we stop to consider data

            for ch in dataVD.keys():
                left = np.searchsorted(
                    dataD[ch][1], timeInit
                )  # index of the first time after or equal to timeInit
                right = np.searchsorted(
                    dataD[ch][1], timeEnd
                )  # index of the time before timeEnd

                dataVD[ch][1] = np.append(dataVD[ch][1], dataD[ch][1][left:right])

        return dataVD


class AtomAnalysis:
    def __init__(self):
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
        ) = "ch0", "ch1", "ch2", "ch3", "ch4", "ch5", "ch6", "ch7"

        self.adt = (
            0.13  # s - minimum atom trapping duration to be considered "good atom"
        )

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

        self.load = DataLoading()

    def dataEv_postSelection(
        self,
        path: str,
        filename: str,
        filetype: str = ".h5",
        kcCounts: int = 2000,
        no=10,
    ):
        print("Post Selecting " + filename + filetype)

        cooling = 25e-3
        photonGate = [0, cooling]

        # ------ We get the data ------
        file = path + "/" + filename + filetype
        filedata = h5py.File(file, mode="r")
        atomnumber = int(len(filedata) / 8)
        pathSave = path + "/goodAtomSelectorFiles/"

        atomList = list(range(0, atomnumber))
        print(atomList)
        # atomList.remove(14)
        # We define and initialize variables before entering the atom loop
        atomDF = pd.DataFrame()
        dataPhoton_grouped = []
        dataTime_grouped = []
        dataPhoton_groupedLC = []
        dataTime_groupedLC = []
        atomIn = []
        atomOut = []
        atomInHisto = []
        atomOutHisto = []
        atomsDuration = []
        atomInNog2 = [x for x in atomIn]
        atomInNog2Histo = [x for x in atomInHisto]
        atomOutNog2 = [x for x in atomOut]
        atomOutNog2Histo = [x for x in atomOutHisto]
        print("atomIn", atomIn)
        print("atomOut", atomOut)

        # ------ We enter the atom loop ------

        for i in tqdm(atomList, file=sys.stdout):
            # Data loading
            # print(i)
            (dataDic) = self.load.data_loading(path, filename, i)

            dataPhotonKC = []
            dataTimeKC = []
            # print(i, dataDic[self.syncFast])
            # print(i, dataDic[self.syncSlow])
            fsdelay = 0.7e-6
            # We calculate the number of photons per trial during the cooling, optical pumping and state detection periods for the Short cavity
            for k, sf in enumerate(dataDic[self.syncFast][1][1:-1]):
                sf = sf + fsdelay
                left1, right1 = np.searchsorted(
                    dataDic[self.kcH][1], [sf + photonGate[0], sf + photonGate[1]]
                )
                left2, right2 = np.searchsorted(
                    dataDic[self.kcV][1], [sf + photonGate[0], sf + photonGate[1]]
                )
                dataPhotonKC.append(right1 - left1 + right2 - left2)
                dataTimeKC.append(sf)

            # no is the number of trials bins that are grouped
            current_dataPhoton_grouped = [
                sum(dataPhotonKC[current : current + no]) / no
                for current in range(0, len(dataPhotonKC) - no, no)
            ]
            dataPhoton_grouped = dataPhoton_grouped + current_dataPhoton_grouped
            current_dataTime_grouped = [
                dataTimeKC[current] for current in range(0, len(dataPhotonKC) - no, no)
            ]
            dataTime_grouped = dataTime_grouped + current_dataTime_grouped
            # print(dataPhotonKC)
            # print(dataTimeKC)
            # print(current_dataPhoton_grouped)
            # print(current_dataTime_grouped)

            dataPhotonLC = []
            dataTimeLC = []

            # We calculate the number of photons per trial for the Long Cavity
            for k, sf in enumerate(dataDic[self.syncFast][1][:-1]):
                sf = sf + fsdelay
                left1, right1 = np.searchsorted(
                    dataDic[self.lcH][1], [sf + photonGate[0], sf + photonGate[1]]
                )
                left2, right2 = np.searchsorted(
                    dataDic[self.lcV][1], [sf + photonGate[0], sf + photonGate[1]]
                )
                dataPhotonLC.append(right1 - left1 + right2 - left2)
                dataTimeLC.append(sf)

            # no = 300  # number of trials bins that are grouped
            current_dataPhoton_groupedLC = [
                sum(dataPhotonLC[current : current + no]) / no
                for current in range(0, len(dataPhotonLC) - no, no)
            ]
            dataPhoton_groupedLC = dataPhoton_groupedLC + current_dataPhoton_groupedLC
            current_dataTime_groupedLC = [
                dataTimeLC[current] for current in range(0, len(dataPhotonLC) - no, no)
            ]
            dataTime_groupedLC = dataTime_groupedLC + current_dataTime_groupedLC

            # print("Mean ph number = ", np.mean(current_dataPhoton_grouped))

            wt_kc = 0.001 * kcCounts  # wt = witness threshold
            wt_lc = -1  # wt = witness threshold
            twot = 4 * kcCounts  # twot = two atom threshold
            inAtom = False
            atomIn_index = 0
            atomOut_index = 0
            for n, j in enumerate(current_dataPhoton_grouped[1:]):
                n = n + 1
                if inAtom == False and j >= wt_kc and j >= wt_lc and j <= twot:
                    atomIn_index = n
                    inAtom = True

                if inAtom:
                    atomOut_index = n

                if inAtom == True and (j < wt_kc or j < wt_lc):
                    atomOut_index = n
                    inAtom == False
                    break

                if inAtom == True and (j > twot):
                    atomOut_index = atomIn_index
                    inAtom == False
                    break

            atomInNog2_index = (
                atomIn_index  # as if there would have been no g2 atom selector
            )
            atomOutNog2_index = (
                atomOut_index  # as if there would have been no g2 atom selector
            )

            try:
                atomIn.append(
                    current_dataTime_grouped[atomIn_index]
                    - dataDic[self.syncSlow][1][0]
                )
                atomInHisto.append(current_dataTime_grouped[atomIn_index])

                atomOut.append(
                    current_dataTime_grouped[atomOut_index]
                    - dataDic[self.syncSlow][1][0]
                )
                atomOutHisto.append(current_dataTime_grouped[atomOut_index])
                atomsDuration.append(atomOut[-1] - atomIn[-1])

                atomInNog2.append(
                    current_dataTime_grouped[atomInNog2_index]
                    - dataDic[self.syncSlow][1][0]
                )
                atomInNog2Histo.append(current_dataTime_grouped[atomInNog2_index])

                atomOutNog2.append(
                    current_dataTime_grouped[atomOutNog2_index]
                    - dataDic[self.syncSlow][1][0]
                )
                atomOutNog2Histo.append(current_dataTime_grouped[atomOutNog2_index])
            except:
                print("Entered in except")
                atomIn.append(0)
                atomInHisto.append(dataDic[self.syncSlow][1][0])
                atomOut.append(0)
                atomOutHisto.append(dataDic[self.syncSlow][1][0])
                atomsDuration.append(atomOut[-1] - atomIn[-1])
                atomInNog2.append(0)
                atomInNog2Histo.append(0)
                atomOutNog2.append(0)
                atomOutNog2Histo.append(0)

        # %% - DATA ALLOCATION IN A DATA FRAME

        # We add the relevant parameters to a data frame
        atomDF["atomsDuration"] = atomsDuration
        atomDF["atomsIn"] = atomIn
        atomDF["atomsOut"] = atomOut
        print("atomIn", atomIn)
        print("atomOut", atomOut)
        # Good atoms are selected, added in the data frame and in a dictionary
        goodAtomsDF = atomDF[(atomDF["atomsDuration"] >= self.adt)]
        goodAtomsDic = {
            i: [goodAtomsDF["atomsIn"][i], goodAtomsDF["atomsOut"][i]]
            for i in list(goodAtomsDF.index)
        }

        # The conditions for good atoms selection are saved in a data frame
        condsDF = pd.DataFrame()
        condsDF["Conditions"] = ["Single atom time threshold (s)"]
        condsDF["Bounds"] = [self.adt]
        atomDicDF = pd.DataFrame.from_dict(
            goodAtomsDic
        )  # data fram with the good atoms dictionary

        # %% ------ We plot the data ------
        plt.close("all")

        f = plt.figure("goodAtomSelector - " + filename, figsize=[17, 14])
        f.suptitle("%s, atom %d, binning = %d" % (filename, i, no))

        ax1 = f.add_subplot(211)
        ax2 = f.add_subplot(212)

        # kc counts plot
        ax1.plot(
            dataTime_grouped,
            dataPhoton_grouped,
            color="tab:orange",
            label="Short Cavity counts",
            ls="None",
            marker=".",
        )
        ax1.vlines(
            atomInHisto, -20, 0, color="grey", linestyle="--", label="atom start time"
        )
        ax1.vlines(
            atomOutHisto, -20, 0, color="red", linestyle="--", label="atom out time"
        )
        ax1.hlines(
            [wt_kc], atomInHisto[0], atomOutHisto[-1], color="tab:green", alpha=0.2
        )
        ax1.hlines([twot], atomInHisto[0], atomOutHisto[-1], color="tab:red", alpha=0.2)
        # lc counts plot
        ax2.plot(
            dataTime_groupedLC,
            dataPhoton_groupedLC,
            color="blue",
            label="Long Cavity counts",
            ls="None",
            marker=".",
        )
        ax2.vlines(
            atomInHisto, -20, 0, color="grey", linestyle="--", label="atom start time"
        )
        ax2.vlines(
            atomOutHisto, -20, 0, color="red", linestyle="--", label="atom out time"
        )
        ax2.hlines(
            [wt_lc], atomInHisto[0], atomOutHisto[-1], color="tab:green", alpha=0.2
        )

        for i in range(len(atomInHisto)):
            if atomOutHisto[i] - atomInHisto[i] >= self.adt:
                ax1.axvspan(
                    atomInHisto[i], atomOutHisto[i], alpha=0.5, color="tab:green"
                )
                ax2.axvspan(
                    atomInHisto[i], atomOutHisto[i], alpha=0.5, color="tab:green"
                )

            if (atomOutHisto[i] == atomInHisto[i]) and (
                atomOutNog2Histo[i] - atomInNog2Histo[i] >= self.adt
            ):
                ax1.axvspan(
                    atomInNog2Histo[i], atomOutNog2Histo[i], alpha=0.3, color="tab:red"
                )
                ax2.axvspan(
                    atomInNog2Histo[i], atomOutNog2Histo[i], alpha=0.3, color="tab:red"
                )

            # print number of atom below
            ax1.text(atomInHisto[i], -20 + 20 * (i % 2), str(i), fontsize=10)
            ax2.text(atomInHisto[i], -20 + 20 * (i % 2), str(i), fontsize=10)

        ax1.set_xlim(xmin=atomInHisto[0] - 2, xmax=atomOutHisto[-1] + 2)
        ax1.legend()

        ax2.set_xlim(xmin=atomInHisto[0] - 2, xmax=atomOutHisto[-1] + 2)
        ax2.legend()

        plt.tight_layout()

        # %% - DATA SAVING
        if self.psSave == True:
            # Check whether the specified path exists or not
            isExist = os.path.exists(pathSave)
            if not isExist:
                os.makedirs(pathSave)

            f.savefig(pathSave + filename + ".png")

            # We save the data frame to an excelfile
            writer = pd.ExcelWriter(pathSave + filename + "_atomParameters.xlsx")
            atomDF.to_excel(writer, sheet_name="atomParameters")
            goodAtomsDF.to_excel(writer, sheet_name="goodAtoms")
            condsDF.to_excel(writer, sheet_name="goodAtomsConds")
            atomDicDF.to_excel(writer, sheet_name="goodAtomsDic")
            writer._save()

            # print(goodAtomsDic)

            # We save the good atoms dictionary
            a_file = open(pathSave + filename + "_goodAtoms.pkl", "wb")
            print(goodAtomsDic)

            list_trappingDuration = []
            for key in goodAtomsDic:
                list_trappingDuration.append(
                    goodAtomsDic[key][1] - goodAtomsDic[key][0]
                )

            print(
                "Average single atom trapping time: (%.2f +/- %.2f)s"
                % (
                    np.mean(list_trappingDuration),
                    np.std(list_trappingDuration)
                    / np.sqrt(np.size(list_trappingDuration)),
                )
            )
            print(
                "Atom trapping probability : %d %%"
                % (len(list_trappingDuration) / len(atomInHisto) * 100)
            )
            print("atom loading trials %.0f" % len(atomInHisto))
            print("Number trapped atoms %.0f" % len(list_trappingDuration))
            print(
                "Duty cycle: %d %%"
                % (
                    sum(list_trappingDuration)
                    / (atomOutHisto[-1] - atomInHisto[0])
                    * 100
                )
            )
            print("Sum trapping duration %.1f" % sum(list_trappingDuration))
            print(
                "Measurement start/end: %.1f/%.1f" % (atomInHisto[0], atomOutHisto[-1])
            )

            pickle.dump(goodAtomsDic, a_file)

            a_file.close()

    def channels_histo(
        self,
        dataDic: dict,
        detectors=["ch2", "ch3", "ch4", "ch7"],
        gates=[0],
        binNum: int = 10000,
        trigger: str = "ch5",
        maxTrigDiff=100e-3,
        fsdelay=[0, 0, 0, 0],
        filename="",
        colors=["grey" for i in range(4)],
    ):
        # dataDic: dictionary with the timestamps generated by data_loading
        # detectors: list of detector strings (e.g. detectors=['lcH','lcV','kcPi','kcV'])
        # gates: list of time gates we want to plot
        # binNum: number of bins
        # trigger: channel used as a trigger
        # maxTrigDiff: maximum time difference between triggers

        syncFast = dataDic[trigger][1]
        histoDic = {
            i: [] for i in detectors
        }  # dictionary where histograms will be stored

        diffFS = np.diff(syncFast)
        print(diffFS)
        # maxFSdur = 11e-3
        maxFSdur = np.amax(
            diffFS[diffFS < maxTrigDiff]
        )  # time difference between atoms are excluded
        print(maxFSdur)
        histotime = np.linspace(0, maxFSdur, binNum)
        binsize = maxTrigDiff / binNum

        for k, det in enumerate(detectors):
            histoDic[det] = np.copy(dataDic[det][1])

            for i in range(len(syncFast) - 1):
                start = syncFast[i] + fsdelay[det]
                FSdur = syncFast[i + 1] - start + fsdelay[det]
                left = np.searchsorted(histoDic[det], start)
                right = np.searchsorted(histoDic[det], start + FSdur)
                histoDic[det][left:right] = histoDic[det][left:right] - start

            histoDic[det] = np.histogram(
                histoDic[det], bins=binNum, range=(0, maxFSdur)
            )[0]

        # Plotting of the histograms
        gateColors = plt.cm.jet(np.linspace(0, 1, len(gates)))  # color map is created

        plt.close(filename + " - Trace Histogram")
        afs = 15

        f = plt.figure(filename + " - Trace Histogram", figsize=[10, 7])

        for i in range(len(detectors)):
            ax = f.add_subplot(len(detectors), 1, i + 1)
            ax.plot(
                histotime * 1e3,
                histoDic[detectors[i]],
                "-",
                label=detectors[i],
                color=colors[i],
            )
            ax1 = ax.twinx()
            ax1.plot(
                histotime * 1e3,
                histoDic[detectors[i]] / (len(syncFast) * binsize),
                "-",
                label=detectors[i],
                color=colors[i],
            )
            ax.legend(loc=6, fontsize=afs)
            ax.set_ylabel("# clicks", fontsize=afs)
            ax1.set_ylabel("rate", fontsize=afs)

            ax.tick_params(axis="y", labelsize=afs)
            ax.tick_params(axis="x", labelsize=afs)

            ax1.tick_params(axis="y", labelsize=afs)

            if i == 0:
                ax1.axhline(
                    y=2400,
                    xmin=histotime[0] * 1e3,
                    xmax=histotime[-1] * 1e3,
                    color="k",
                    ls="--",
                )  # 2-2' pumping dark counts

            if i == 1:
                ax1.axhline(
                    y=1600,
                    xmin=histotime[0] * 1e3,
                    xmax=histotime[-1] * 1e3,
                    color="k",
                    ls="--",
                )  # 2-2' pumping dark counts

            for i, gate in enumerate(gates):
                ax.axvline(gate * 1e3, color=gateColors[i])
        ax.set_xlabel("Time (ms)", fontsize=afs)

        f.suptitle(filename + " - Trace Histogram")

        plt.tight_layout()

        return f

    def dataEval_noramlModeSpectroscopy(
        self,
        path: str,
        file_list: List[str],
        filetype: str = "h5",
        ParamDic={},
    ):
        """
        INPUTS
        path: path to directory where to find the dataDic
        filename: name of the file to be analyzed
        filetype: h5 default

        OUTPUT:
        Whatever
        """
        print("\nAnalysing: ")

        for filename in file_list:
            print(filename + ", ")

            file_postSelected = (
                path + "/goodAtomSelectorFiles/" + filename + "_goodAtoms.pkl"
            )
            if not os.path.exists(file_postSelected):
                self.dataEv_postSelection(
                    path, filename, filetype, kcCounts=2000, no=10
                )

            a_file = open(file_postSelected, "rb")
            atomDic = pickle.load(a_file)

            print(atomDic)

            # ------ We get the data ------
            dataDic = self.load.data_goodAtoms(path, filename, atomDic)

            # sequence parameters
            trigger_delay = 3.15e-6
            coolDur = 400e-6
            OPdur = 200e-6
            pulse_delay = 33.5e-6
            pulseDur = 7e-6

            # sequence gates
            opgate = [trigger_delay + coolDur, trigger_delay + coolDur + OPdur]
            writegate = [opgate[1] + pulse_delay, opgate[1] + pulse_delay + pulseDur]
            gates = writegate

            # plot histogram
            seqDur = 0.7e-3
            print(seqDur)
            trig = self.syncFast
            scanTrig = self.syncFast2
            maxTrigDiff = seqDur  # should be bigger than triggers time difference
            binsize = 20 * 1e-9
            binNum = int(maxTrigDiff / binsize)
            detectors = [self.kcH, self.kcV, self.lcH, self.lcV, self.sdTrig]
            colors = ["violet", "violet", "tab:blue", "tab:blue", "orange"]
            fsdelay = {
                self.kcH: 0,
                self.kcV: 12e-9,
                self.lcH: 0,
                self.lcV: 0.0,
                self.sdTrig: 0.0,
            }
            # chfig = self.channels_histo(
            #     dataDic,
            #     detectors,
            #     gates,
            #     binNum,
            #     trig,
            #     maxTrigDiff,
            #     fsdelay,
            #     filename,
            #     colors,
            # )
            # plt.show(block=True)
            #
            # extract data for Normal Mode Spectrum
            fsDiff = dataDic[trig][1][-1] - dataDic[trig][1][-2]
            scan_duration = (
                fsDiff * ParamDic["PointsPerScan"] * ParamDic["TrialsPerPoint"]
            )
            print("Fast trigger period = %.5f s" % (scan_duration))
            print("Scan period = %.5f s" % (scan_duration))
            freq_NMS = (
                np.linspace(0, ParamDic["freqSpan"], int(ParamDic["PointsPerScan"] / 2))
                - ParamDic["freqCenter"]
            )
            R_binary_up = dict(
                (i, []) for i in range(int(ParamDic["PointsPerScan"] / 2))
            )
            R_binary_down = dict(
                (i, []) for i in range(int(ParamDic["PointsPerScan"] / 2))
            )
            T_binary_up = dict(
                (i, []) for i in range(int(ParamDic["PointsPerScan"] / 2))
            )
            T_binary_down = dict(
                (i, []) for i in range(int(ParamDic["PointsPerScan"] / 2))
            )
            NoPoints = int(ParamDic["PointsPerScan"] / 2)
            NoTrials = ParamDic["TrialsPerPoint"]
            print("""
                =================
                Analyse atom data
                =================
                --- Looping over Fast Sequence triggers ---     
                """)
            # Loop over all scan trigger
            for i in range(len(dataDic[scanTrig][1][:-1])):
                # remove triggers of incomplete scan
                sfsDiff = dataDic[scanTrig][1][i + 1] - dataDic[scanTrig][1][i]
                if sfsDiff > scan_duration * 1.1:
                    print(sfsDiff)
                    np.delete(dataDic[scanTrig][1], i)
                else:
                    # start of each spectrum
                    indexSF0 = np.searchsorted(
                        dataDic[trig][1], dataDic[scanTrig][1][i]
                    )
                    for point in range(0, int(ParamDic["PointsPerScan"] / 2), 1):
                        # iterate through all trials per step
                        for trial in range(0, ParamDic["TrialsPerPoint"], 1):
                            sf = indexSF0 + point * ParamDic["TrialsPerPoint"] + trial
                            # Herald Transmission
                            left_H, right_H = np.searchsorted(
                                dataDic[self.kcH][1],
                                [
                                    dataDic[trig][1][sf] + writegate[0],
                                    dataDic[trig][1][sf] + writegate[1],
                                ],
                            )
                            left_V, right_V = np.searchsorted(
                                dataDic[self.kcV][1],
                                [
                                    dataDic[trig][1][sf] + writegate[0],
                                    dataDic[trig][1][sf] + writegate[1],
                                ],
                            )
                            T_binary_up[point].append(
                                right_H - left_H + right_V - left_V
                            )

                            # Herald Transmission
                            left_H, right_H = np.searchsorted(
                                dataDic[self.kcH][1],
                                [
                                    dataDic[trig][1][sf + NoPoints * NoTrials]
                                    + writegate[0],
                                    dataDic[trig][1][sf + NoPoints * NoTrials]
                                    + writegate[1],
                                ],
                            )
                            left_V, right_V = np.searchsorted(
                                dataDic[self.kcV][1],
                                [
                                    dataDic[trig][1][sf + NoPoints * NoTrials]
                                    + writegate[0],
                                    dataDic[trig][1][sf + NoPoints * NoTrials]
                                    + writegate[1],
                                ],
                            )
                            T_binary_down[point].append(
                                right_H - left_H + right_V - left_V
                            )

                            # Reflection
                            left_H, right_H = np.searchsorted(
                                dataDic[self.lcH][1],
                                [
                                    dataDic[trig][1][sf] + writegate[0],
                                    dataDic[trig][1][sf] + writegate[1],
                                ],
                            )
                            left_V, right_V = np.searchsorted(
                                dataDic[self.lcV][1],
                                [
                                    dataDic[trig][1][sf] + writegate[0],
                                    dataDic[trig][1][sf] + writegate[1],
                                ],
                            )
                            R_binary_up[point].append(
                                right_V - left_V + right_H - left_H
                            )

                            left_H, right_H = np.searchsorted(
                                dataDic[self.lcH][1],
                                [
                                    dataDic[trig][1][sf + 60 * 20] + writegate[0],
                                    dataDic[trig][1][sf + 60 * 20] + writegate[1],
                                ],
                            )
                            left_V, right_V = np.searchsorted(
                                dataDic[self.lcV][1],
                                [
                                    dataDic[trig][1][sf + 60 * 20] + writegate[0],
                                    dataDic[trig][1][sf + 60 * 20] + writegate[1],
                                ],
                            )
                            R_binary_down[point].append(
                                right_V - left_V + right_H - left_H
                            )

        list_SDmean = []
        list_err_SDmean = []
        for i in range(2):
            list_SDmean.append(np.zeros(int(ParamDic["PointsPerScan"] / 2)))
            list_err_SDmean.append(np.zeros(int(ParamDic["PointsPerScan"] / 2)))

        for point in range(0, int(ParamDic["PointsPerScan"] / 2), 1):
            list_SDmean[0][point] = np.mean(T_binary_up[point])
            list_err_SDmean[0][point] = np.std(T_binary_up[point]) / np.sqrt(
                len(T_binary_up[point])
            )

            list_SDmean[1][point] = np.mean(T_binary_down[point])
            list_err_SDmean[1][point] = np.std(T_binary_down[point]) / np.sqrt(
                len(T_binary_down[point])
            )

        list_SDmean[0] = list_SDmean[0] + np.flip(list_SDmean[1])
        list_err_SDmean[0] = list_err_SDmean[0] + np.flip(list_err_SDmean[1])

        # ------ Plot for Photon ----
        # 1. Define the exponentially modified Gaussian (EMG)
        def Gaussian(x, A=1, mu=0, sigma=1):
            """Gaussian normal distribution using NumPy."""
            return (A * 1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(
                -((x - mu) ** 2) / (2 * sigma**2)
            )

        def R_coupled(detuning, f_res, g, c):
            """Reflection model of a coupled atom–cavity system.

            Parameters
            ----------
            detuning : float or array
                Probe–cavity detuning (independent variable, e.g. MHz).
            f_res : float
                Effective resonance frequency offset of the coupled system.
            A : float
                Overall amplitude / normalization factor.
            g : float
                Atom–cavity coupling rate.
            kappa : float
                Total cavity field decay rate (linewidth).
            kappa_oc : float
                Cavity outcoupling rate (through monitored mirror/fiber).
            MM_rf : float
                Reflection–fiber mode matching efficiency.
            MM_fc : float
                Fiber–cavity mode matching efficiency.
            gamma : float
                Free-space atomic decay rate.
            offset : float
                Constant vertical offset (background level).
            a : float
                Linear slope correction, making the effective decay rate
                detuning-dependent via Gamma = gamma + a*detuning.

            Returns
            -------
            float or array
                Reflected intensity at the given detuning.
            """
            # Effective atomic decay including optional detuning-dependent slope
            # Gamma = gamma + a * detuning
            # Gamma = gamma

            # Cooperativity-like denominator term
            C_d = g**2 / (
                2 * (3.0 + 1j * (detuning - f_res)) * (58 + 1j * (detuning - f_res))
            )

            # Reflection amplitude → squared magnitude → intensity
            return (
                abs(
                    0.978
                    - (0.871 * np.exp(1j * 0.0) ** 2)
                    * 2
                    * 58
                    * 0.85
                    / ((58 + 1j * (detuning - f_res)) * (2 * C_d + 1))
                )
                ** 2
                + c
            )

        # 2. Fit the EMG to the histogram

        p0 = [
            -40,  # f_res (MHz offset of resonance)
            50,  # g (coupling strength, MHz)
            0,
        ]

        bounds = (
            [-100, 10, -1],
            [0, 50, 1],
        )  # Avoid zero/negative sigma/lambda

        norm_coefficient = 0.20780113085143606

        popt, pcov = curve_fit(
            R_coupled,
            freq_NMS,
            list_SDmean[0] / norm_coefficient,
            p0=p0,
            bounds=bounds,
            maxfev=100000,
        )

        # print(popt)
        pcov = np.sqrt(np.diag(pcov))
        plt.rcParams.update({"font.size": 14})
        phfig = plt.figure(figsize=[12, 8])
        phfig.suptitle(
            # filename + "\n Memory Spectroscopy with KC @ +500MHz detuning"
            # "Normal Mode Spectroscopy with KC\n"
            r"g: %.1f MHz $\pm$ %.1f MHz" % (popt[1], pcov[1])
        )
        ax = phfig.add_subplot(1, 1, 1)
        ax.errorbar(
            freq_NMS,
            list_SDmean[0] / norm_coefficient,
            list_err_SDmean[0] / norm_coefficient,
            linestyle="",
            marker="o",
            label="Measurement data",
        )
        ax.plot(
            freq_NMS,
            R_coupled(freq_NMS, *popt),
            label="Model fit",
            color="red",
            linewidth=3,
            linestyle="-.",
        )
        # ax.plot(
        #     freq_NMS,
        #     R_coupled(freq_NMS, *p0),
        #     label="Personal guess",
        #     color="green",
        #     linewidth=3,
        #     linestyle="-.",
        # )
        ax.set_ylabel("Reflection (a.u.)")
        ax.legend()
        ax.set_xlabel(r"Detuning $\Delta$ (MHz)")
        plt.tight_layout()
        plt.show()
        phfig.savefig(f"{path}/{filename}_reflection_spectrum.svg")
        # with open("data_from_fit.csv", "w") as file:
        #     file.write("x,y\n")
        #     for i in range(len(freq_NMS)):
        #         file.write(f"{freq_NMS[i]},{R_coupled(freq_NMS[i], *popt)}\n")


if __name__ == "__main__":
    """
    This script analyzes normal mode spectroscopy measurement data
    """

    # ------ Setting the print options for numpy ------
    np.set_printoptions(legacy="1.25")

    # ------ Setting flags------
    loadFromDir: bool = False  # if True, counts dictionary is loaded from the directory

    # ------ Definition of data sources and destinations------
    path: str = os.path.abspath(os.path.dirname(__file__))
    file_list = ["24_09_25_KC_Spectroscopy_2_1_pi_75_375_MHz_100_points_3"]
    filetype: str = ".h5"

    # ------ Printout of flags and measurement ------
    print("Load coincidences from pickle file in Directory : ", loadFromDir)

    # ------ Begin Analysis ------
    analysis = AtomAnalysis()

    ParamDic = {
        "freqSpan": 250,  # in MHz
        "PointsPerScan": 200,  # including up and down ramp
        "TrialsPerPoint": 40,
        "freqCenter": 125,
    }

    analysis.dataEval_noramlModeSpectroscopy(
        path, file_list, filetype, ParamDic=ParamDic
    )
