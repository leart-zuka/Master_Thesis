import os
import sys
import time
import h5py
from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib import colormaps
import numpy as np
import lmfit
from scipy.optimize import curve_fit, minimize
import pandas as pd
import pickle
from openpyxl import load_workbook
from datetime import datetime
import matplotlib.colors as mcolors
from typing import List
from qutip import *
from qutip import gates


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

    def atomDicCleaner(self, atomDic: dict, specDuration):
        for key in atomDic:
            specNr = (atomDic[key][1] - atomDic[key][0]) // (specDuration)
            atomDic[key][1] = specNr * specDuration + atomDic[key][0]
        return atomDic


# --- Functions for curve fitting ---
class fittingFunctions:
    def __init__(self):
        print("fittingFunctions loaded")
        return None

    def exponential(self, x, A0, tau, offset):
        return A0 * np.exp(-x / tau) + offset

    def rabiRotation(self, time, freq, decRate, A, off):
        """
        INPUTS
        time: duration of the transfer pulse
        rabiFreq: Rabi frequency
        decRate:  oscillations amplitude exponential decay rate
        OUTPUTS
        p1: probability for the atom to be in state 1 (F=2)
        """
        # p1 = (A * np.exp(-decRate * 1e3 * time) * (np.sin(2 * np.pi * freq * 1e3 * time)) ** 2 + A / 2 * (
        #             1 - np.exp(-decRate * 1e3 * time))) + off

        p1 = (
            A
            * (
                (np.sin(2 * np.pi * freq * 1e3 * time) ** 2 - 0.5)
                * np.exp(-decRate * 1e3 * time)
                + 0.5
            )
            + off
        )

        return p1

    def ramseyRotation(self, time, freq, decRate, A, off):
        """
        INPUTS
        time: duration of the transfer pulse
        rabiFreq: Rabi frequency
        decRate:  oscillations amplitude exponential decay rate
        OUTPUTS
        p1: probability for the atom to be in state 1 (F=2)
        """

        p1 = (
            A
            * (
                (np.cos(2 * np.pi * freq * 1e3 * time) ** 2 - 0.5)
                * np.exp(-decRate * 1e3 * time)
                + 0.5
            )
            + off
        )

        return p1

    def lorenzian(self, freq, amp1, foff, k1, offset):

        p = amp1 * k1**2 / ((freq - foff) ** 2 + k1**2) + offset

        return p


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
        no=1,
    ):
        print("Post Selecting " + filename + filetype)

        cooling = 25e-3
        photonGate = [0, cooling]

        # ------ We get the data ------
        file = path + "/" + filename + filetype
        filedata = h5py.File(file, mode="r")
        atomnumber = int(len(filedata) / 8)
        pathSave = path + "\\goodAtomSelectorFiles\\"

        atomList = range(0, atomnumber)

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

        # ------ We enter the atom loop ------

        for i in tqdm(atomList, file=sys.stdout):
            # Data loading

            # (dataDic) = gf.data_loading(path, filename, i)
            (dataDic) = self.load.data_loading(path, filename, i)

            dataPhotonKC = []
            dataTimeKC = []

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

            # no = 300  # number of trials bins that are grouped
            current_dataPhoton_grouped = [
                sum(dataPhotonKC[current : current + no]) / no
                for current in range(0, len(dataPhotonKC) - no, no)
            ]
            dataPhoton_grouped = dataPhoton_grouped + current_dataPhoton_grouped
            current_dataTime_grouped = [
                dataTimeKC[current] for current in range(0, len(dataPhotonKC) - no, no)
            ]
            dataTime_grouped = dataTime_grouped + current_dataTime_grouped

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

            wt_kc = 0.3 * kcCounts  # wt = witness threshold
            wt_lc = -1  # wt = witness threshold
            twot = 1.4 * kcCounts  # twot = two atom threshold
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
            print(
                "Duty cycle: %d %%"
                % (
                    sum(list_trappingDuration)
                    / (atomOutHisto[-1] - atomInHisto[0])
                    * 100
                )
            )

            pickle.dump(goodAtomsDic, a_file)

            a_file.close()

    def dataEv_atomPhotonCoincidences(
        self,
        path: str,
        file_list: List[str],
        filetype: str = ".h5",
        position: int = 0,
        heralded: bool = True,
    ):
        # ------ Definition of data sources ------

        """analysing files"""

        # ------ Definition of parameters ------

        pathSave = path + "\\atomPhotonCoincidenceFiles\\"
        efficiencyDF = pd.DataFrame()
        sdDF = pd.DataFrame()
        g2DF = pd.DataFrame()

        tlcH, tlcV, tkcH, tkcV = [], [], [], []
        photonlcH, photonlcV, photonkcH, photonkcV = [], [], [], []
        photonSD = [[], []]  # n of photons in SD gate for each trial

        kclcHF1, kclcHF2, kclcVF1, kclcVF2 = 0, 0, 0, 0
        lcHF1, lcHF2, lcVF1, lcVF2 = 0, 0, 0, 0
        fstot = 0
        sdTrigtot = 0

        print("\nAnalysing: ")
        for filename in file_list:
            print(filename + ", ")
        print("trigger ", position + 1)

        for filename in file_list:
            file_postSelected = (
                path + "\\goodAtomSelectorFiles\\" + filename + "_goodAtoms.pkl"
            )
            if not os.path.exists(file_postSelected):
                self.dataEv_postSelection(path, filename, filetype)

            a_file = open(file_postSelected, "rb")
            atomList = pickle.load(a_file)

            print(atomList)

            # ------ We get the data ------
            # dataDic = gf.data_goodAtoms(path, filename, atomList)
            dataDic = self.load.data_goodAtoms(path, filename, atomList)

            cooldur = 20.006551e-3 + 20e-3 - 2e-6 + 1.170e-6
            OPdur = 150e-6
            gap = 1 * 1e-6
            excdur = 5e-6
            MW = 14.12e-6 + 3e-6 + 21.8e-6
            SD = 7.5e-6

            wait = 1000 * 1e-6

            coolgate = [0, cooldur]
            opgate = [cooldur, cooldur + OPdur]

            excgate = [
                opgate[1] + gap + 7.38 * 1e-6,
                opgate[1] + gap + 7.38 * 1e-6 + 600e-9,
            ]
            MWgate = [
                opgate[1] + 3 * gap + excdur,
                opgate[1] + 3 * gap + excdur + MW + 50e-6,
            ]

            # SD and coincidence gates after "wait"  time
            SDgate = [40.26886e-3, 40.26886e-3 + SD]
            coincgate = [40.265e-3, 40.27e-3]
            hgate = [-20e-9, +20e-9]

            gates = excgate + SDgate + coincgate
            seqDur = 40e-3

            # trigger = self.syncFast
            # maxTrigDiff = seqDur + 10e-6 + 50e-6# should be bigger than triggers time difference
            trigger = self.syncFast
            maxTrigDiff = (
                seqDur + 20e-3
            )  # should be bigger than triggers time difference
            binsize = 0.5 * 1e-6
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
            trfig = self.channels_histo(
                dataDic,
                detectors,
                gates,
                binNum,
                trigger,
                maxTrigDiff,
                fsdelay,
                filename,
                colors,
            )
            SDthreshold = (
                1  # number of photons starting from which a bright state is witnessed
            )

            print("""
            ============
            Atoms data
            ============
                 
            --- Looping over Fast Sequence triggers ---     
            """)

            # trig = self.syncFast
            trig = self.syncFast2

            time.sleep(0.1)
            # Loop over all the fast sequence triggers
            for fsr in tqdm(dataDic[trig][1][:-5], file=sys.stdout):
                fs = dataDic[self.syncFast][1][
                    np.searchsorted(dataDic[self.syncFast][1], fsr) + position
                ]  # + position*1e-6

                # detector lcH
                sf = fs + fsdelay[self.lcH]
                left, right = np.searchsorted(
                    dataDic[self.lcH][1], [sf + excgate[0], sf + excgate[1]]
                )
                (
                    tlcH.extend(dataDic[self.lcH][1][left:right] - sf)
                    if left - right != 0
                    else None
                )
                photonlcH.append(right - left)
                # detector lcV
                sf = fs + fsdelay[self.lcV]
                left, right = np.searchsorted(
                    dataDic[self.lcV][1], [sf + excgate[0], sf + excgate[1]]
                )
                (
                    tlcV.extend(dataDic[self.lcV][1][left:right] - sf)
                    if left - right != 0
                    else None
                )
                photonlcV.append(right - left)
                # kcH herald
                sf = fs + fsdelay[self.kcH]
                left, right = np.searchsorted(
                    dataDic[self.kcH][1], [sf + excgate[0], sf + excgate[1]]
                )
                (
                    tkcH.extend(dataDic[self.kcH][1][left:right] - sf)
                    if left - right != 0
                    else None
                )
                photonkcH.append(right - left)
                # kcV herald
                sf = fs + fsdelay[self.kcV]
                left, right = np.searchsorted(
                    dataDic[self.kcV][1], [sf + excgate[0], sf + excgate[1]]
                )
                (
                    tkcV.extend(dataDic[self.kcV][1][left:right] - sf)
                    if left - right != 0
                    else None
                )
                photonkcV.append(right - left)

                # # --- using only 20ns time gate depending on Herald detection ---
                # # kcH herald
                # sf = fs + fsdelay[self.kcH]
                # left, right = np.searchsorted(dataDic[self.kcH][1], [sf + excgate[0], sf + excgate[1]])
                # if right - left != 0:
                #     tH = dataDic[self.kcH][1][left:right] - sf
                #     tkcH.extend(tH)
                # else:
                #     tH = None
                # photonkcH.append(right - left)
                # # kcV herald
                # sf = fs + fsdelay[self.kcV]
                # left, right = np.searchsorted(dataDic[self.kcV][1], [sf + excgate[0], sf + excgate[1]])
                # if right - left != 0:
                #     tV = dataDic[self.kcV][1][left:right] - sf
                #     tkcV.extend(tV)
                # else:
                #     tV = None
                # photonkcV.append(right - left)

                # if tV is not None:
                #     # detector lcH
                #     sf = fs + tV[0] + fsdelay[self.lcH]
                #     left, right = np.searchsorted(dataDic[self.lcH][1], [sf + hgate[0], sf + hgate[1]])
                #     (tlcH.extend(dataDic[self.lcH][1][left:right] - sf) if left - right != 0 else None)
                #     photonlcH.append(right - left)
                #     # detector lcV
                #     sf = fs + tV[0] + fsdelay[self.lcV]
                #     left, right = np.searchsorted(dataDic[self.lcV][1], [sf + hgate[0], sf + hgate[1]])
                #     (tlcV.extend(dataDic[self.lcV][1][left:right] - sf) if left - right != 0 else None)
                #     photonlcV.append(right - left)
                # if tH is not None:
                #     # detector lcH
                #     sf = fs + tH[0] + fsdelay[self.lcH]
                #     left, right = np.searchsorted(dataDic[self.lcH][1], [sf + hgate[0], sf + hgate[1]])
                #     (tlcH.extend(dataDic[self.lcH][1][left:right] - sf) if left - right != 0 else None)
                #     photonlcH.append(right - left)
                #     # detector lcV
                #     sf = fs + tH[0] + fsdelay[self.lcV]
                #     left, right = np.searchsorted(dataDic[self.lcV][1], [sf + hgate[0], sf + hgate[1]])
                #     (tlcV.extend(dataDic[self.lcV][1][left:right] - sf) if left - right != 0 else None)
                #     photonlcV.append(right - left)

                # if tV is None and tH is None:
                #     # detector lcH
                #     sf = fs + fsdelay[self.lcH]
                #     left, right = np.searchsorted(dataDic[self.lcH][1], [sf + excgate[0], sf + excgate[1]])
                #     photonlcH.append(right - left)
                #     # detector lcV
                #     sf = fs + fsdelay[self.lcV]
                #     left, right = np.searchsorted(dataDic[self.lcV][1], [sf + excgate[0], sf + excgate[1]])
                #     photonlcV.append(right - left)
                # elif tV is not None:
                #     # detector lcH
                #     sf = fs + tV[0] + fsdelay[self.lcH]
                #     left, right = np.searchsorted(dataDic[self.lcH][1], [sf + hgate[0], sf + hgate[1]])
                #     (tlcH.extend(dataDic[self.lcH][1][left:right] - sf) if left - right != 0 else None)
                #     photonlcH.append(right - left)
                #     # detector lcV
                #     sf = fs + tV[0] + fsdelay[self.lcV]
                #     left, right = np.searchsorted(dataDic[self.lcV][1], [sf + hgate[0], sf + hgate[1]])
                #     (tlcV.extend(dataDic[self.lcV][1][left:right] - sf) if left - right != 0 else None)
                #     photonlcV.append(right - left)
                # elif tH is not None:
                #     # detector lcH
                #     sf = fs + tH[0] + fsdelay[self.lcH]
                #     left, right = np.searchsorted(dataDic[self.lcH][1], [sf + hgate[0], sf + hgate[1]])
                #     (tlcH.extend(dataDic[self.lcH][1][left:right] - sf) if left - right != 0 else None)
                #     photonlcH.append(right - left)
                #     # detector lcV
                #     sf = fs + tH[0] + fsdelay[self.lcV]
                #     left, right = np.searchsorted(dataDic[self.lcV][1], [sf + hgate[0], sf + hgate[1]])
                #     (tlcV.extend(dataDic[self.lcV][1][left:right] - sf) if left - right != 0 else None)
                #     photonlcV.append(right - left)

                # SD
                sf = fs
                left_h, right_h = np.searchsorted(
                    dataDic[self.kcH][1], [sf + SDgate[0], sf + SDgate[1]]
                )
                left_v, right_v = np.searchsorted(
                    dataDic[self.kcV][1], [sf + SDgate[0], sf + SDgate[1]]
                )
                left_cond, right_cond = np.searchsorted(
                    dataDic[self.sdTrig][1], [sf + coincgate[0], sf + coincgate[1]]
                )

                if right_cond - left_cond >= 1:
                    sdTrigtot += 1

                    photonSD[0].append((right_v - left_v) + (right_h - left_h))
                    (
                        photonSD[1].extend(dataDic[self.kcH][1][left_h:right_h] - sf)
                        if left_h - right_h != 0
                        else None
                    )
                    (
                        photonSD[1].extend(dataDic[self.kcV][1][left_v:right_v] - sf)
                        if left_v - right_v != 0
                        else None
                    )

                    # heralded coincidences
                    # if ((photonkcH[-1] or photonkcV[-1])):
                    if photonkcH[-1]:
                        kclcHF1 += photonlcH[-1] & (photonSD[0][-1] == 0)
                        kclcHF2 += photonlcH[-1] & (photonSD[0][-1] >= SDthreshold)
                        kclcVF1 += photonlcV[-1] & (photonSD[0][-1] == 0)
                        kclcVF2 += photonlcV[-1] & (photonSD[0][-1] >= SDthreshold)

                    # non heralded coincidences
                    if photonlcH[-1] or photonlcV[-1]:
                        lcHF1 += photonlcH[-1] & (photonSD[0][-1] == 0)
                        lcHF2 += photonlcH[-1] & (photonSD[0][-1] >= SDthreshold)
                        lcVF1 += photonlcV[-1] & (photonSD[0][-1] == 0)
                        lcVF2 += photonlcV[-1] & (photonSD[0][-1] >= SDthreshold)

                fstot += 1

        # %% --- Plotting, Fitting and Saving

        print("tot photons = ", (sum(photonlcH) + sum(photonlcV)))
        photonlcEff = (sum(photonlcH) + sum(photonlcV)) / fstot
        photonlcEfferr = np.sqrt(sum(photonlcH) + sum(photonlcV)) / fstot
        photonkcEff = (sum(photonkcH) + sum(photonkcV)) / fstot
        photonkcEfferr = np.sqrt(sum(photonkcH) + sum(photonkcV)) / fstot
        # eff_kc_upon_lc = sum((np.array(photonkcH)+np.array(photonkcV))*(np.array(photonlcH)+np.array(photonlcV)))/(sum(photonlcH) + sum(photonlcV))
        # eff_kc_upon_lcerr = np.sqrt(sum((np.array(photonkcH)+np.array(photonkcV))*(np.array(photonlcH)+np.array(photonlcV))))/(sum(photonlcH) + sum(photonlcV))
        # eff_lc_upon_kc = sum((np.array(photonkcH)+np.array(photonkcV))*(np.array(photonlcH)+np.array(photonlcV)))/(sum(photonkcH) + sum(photonkcV))
        # eff_lc_upon_kcerr = np.sqrt(sum((np.array(photonkcH)+np.array(photonkcV))*(np.array(photonlcH)+np.array(photonlcV))))/(sum(photonkcH) + sum(photonkcV))
        # coinc_n = sum((np.array(photonkcH)+np.array(photonkcV))*(np.array(photonlcH)+np.array(photonlcV)))

        print("\nnumber of atom fs triggers = ", fstot)
        efflcText = "Photon detection effciency lC = %.2f pm %.2f %%" % (
            photonlcEff * 100,
            photonlcEfferr * 100,
        )
        print(efflcText)
        effkcText = "Photon detection effciency kC = %.2f pm %.2f %%" % (
            photonkcEff * 100,
            photonkcEfferr * 100,
        )
        print(effkcText)
        # eff_lc_upon_kc_Text = "Photon detection effciency lC upon kc detection = %.2f pm %.2f%%"%(eff_lc_upon_kc*100, eff_lc_upon_kcerr*100)
        # print(eff_lc_upon_kc_Text)
        # eff_kc_upon_lc_Text = "Photon detection effciency kC upon lc detection = %.2f pm %.2f%%"%(eff_kc_upon_lc*100, eff_kc_upon_lcerr*100)
        # print(eff_kc_upon_lc_Text)
        # print("\n")
        # print("Number of lc-kc coincidences = %d"%(coinc_n))
        print("Number of SD triggers = %d" % (sdTrigtot))
        # print("lc-kc coincidences probability = %.2f pm %.2f%%"%(coinc_n/fstot*100, np.sqrt(coinc_n)/fstot * 100))

        # ------ Plot for Photon ----
        bin_width = 5e-9
        # historange = [excgate[0], excgate[1]]
        historange = [hgate[0], hgate[1]]
        histoBinNum = int((historange[1] - historange[0]) / bin_width)

        histogram_lc = np.histogram(tlcH + tlcV, range=historange, bins=histoBinNum)[0]
        histogram_kc = np.histogram(tkcH + tkcV, range=historange, bins=histoBinNum)[0]

        plt.rcParams.update({"font.size": 20})
        spttl = filename + "_position" + str(position) + " - Photon"
        plt.close(spttl)
        phfig = plt.figure(spttl, figsize=[10, 7])
        phfig.suptitle(
            filename
            + "\nbin width= %.1f ns" % (bin_width * 1e9)
            + "\n"
            + efflcText
            + "\n"
            + effkcText
        )
        ax1 = phfig.add_subplot(2, 1, 1)
        histotime = (
            np.linspace(historange[0], historange[1], histoBinNum) - historange[0]
        ) * 1e9  # ns
        ax1.bar(
            histotime,
            histogram_lc / fstot,
            width=bin_width * 1e9,
            color=mcolors.CSS4_COLORS["cornflowerblue"],
            label="lcPhoton",
        )
        ax = phfig.add_subplot(2, 1, 2)
        ax.bar(
            histotime,
            histogram_kc / fstot,
            width=bin_width * 1e9,
            color="orange",
            label="kcPhoton",
        )

        ax.legend()
        ax.set_xlabel("time (ns)")
        ax.set_ylabel("counts/trials")
        plt.tight_layout()

        ### ------ Plot for SD ----
        bin_width = 0.25 * 1e-6
        gap = 1e-6
        historange = [SDgate[0] - gap, SDgate[1] + gap]
        histoBinNum = int((historange[1] - historange[0]) / bin_width)

        histogram_sd = np.histogram(photonSD[1], range=historange, bins=histoBinNum)[0]

        plt.rcParams.update({"font.size": 20})

        spttl = filename + "_position" + str(position) + " - State Detection"
        plt.close(spttl)
        sdfig = plt.figure(spttl, figsize=[10, 7])
        sdfig.suptitle(
            filename
            + "\nbin width= %.1f ns" % (bin_width * 1e9)
            + "\n"
            + efflcText
            + "\n"
            + effkcText
        )
        ax = sdfig.add_subplot(1, 1, 1)
        histotime = (
            np.linspace(historange[0], historange[1], histoBinNum) - historange[0]
        ) * 1e9  # ns
        ax.bar(
            histotime,
            histogram_sd / fstot,
            width=bin_width * 1e9,
            color="orange",
            label="SD",
        )
        ax.vlines(
            np.array([gap, gap + (SDgate[1] - SDgate[0])]) * 1e9,
            ymin=0,
            ymax=max(histogram_sd / fstot),
            ls="--",
            color="k",
        )

        ax.legend()
        ax.set_xlabel("time (ns)")
        ax.set_ylabel("counts/trials")

        # %% ---- lc photon g2 autocorrelation analysis ----
        trialShifts = np.arange(-10, 10)
        coincidences_kc = np.zeros(len(trialShifts))
        coincidences_lc = np.zeros(len(trialShifts))

        for i, shift in enumerate(trialShifts):
            coincidences_kc[i] = np.inner(photonkcH, np.roll(photonkcV, shift))
            coincidences_lc[i] = np.inner(photonlcH, np.roll(photonlcV, shift))

        g2_kc = coincidences_kc / np.average(np.delete(coincidences_kc, 10))
        g2_lc = coincidences_lc / np.average(np.delete(coincidences_lc, 10))
        errg2_kc = (
            np.clip(coincidences_kc, 1, 1e15)
            / np.average(np.delete(coincidences_kc, 10))
            * np.sqrt(
                1 / np.clip(coincidences_kc, 1, 1e15)
                + 1 / np.sum(np.delete(coincidences_kc, 10))
            )
        )
        errg2_lc = (
            np.clip(coincidences_lc, 1, 1e15)
            / np.average(np.delete(coincidences_lc, 10))
            * np.sqrt(
                1 / np.clip(coincidences_lc, 1, 1e15)
                + 1 / np.sum(np.delete(coincidences_lc, 10))
            )
        )

        print("""
        ============
        g2 Analysis
        ============    
        """)

        print(
            "Photons in kcPi = %d\nPhotons in kcV = %d"
            % (sum(photonkcH), sum(photonkcV))
        )
        print("coincidences(0) kC = ", coincidences_kc[10])
        kcText = "kc: g2(0) = %.2f pm %.2f" % (g2_kc[10], errg2_kc[10])
        print(kcText)
        print(
            "\nPhotons in det1 = %d\nPhotons in det2 = %d"
            % (sum(photonlcH), sum(photonlcV))
        )
        print("coincidences(0) kC = ", coincidences_lc[10])
        lcText = "lc: g2(0) = %.2f pm %.2f" % (g2_lc[10], errg2_lc[10])
        print(lcText)

        spttl = filename + "_position" + str(position) + " - g2"
        plt.close(spttl)
        g2fig = plt.figure(spttl, figsize=[10, 7])
        ax1 = g2fig.add_subplot(2, 1, 1)
        ax2 = g2fig.add_subplot(2, 1, 2)

        g2fig.suptitle(kcText + "\n" + lcText)
        ax2.bar(
            trialShifts * 40,
            g2_kc,
            yerr=errg2_kc,
            width=8,
            color="orange",
            label="kcPhoton",
        )
        ax1.bar(
            trialShifts * 40,
            g2_lc,
            yerr=errg2_lc,
            width=8,
            color=mcolors.CSS4_COLORS["cornflowerblue"],
            label="lcPhoton",
        )

        ax2.set_xlabel(r"$\tau$ (ms)")

        plt.tight_layout()

        if self.dataSave:
            # Check whether the specified path exists or not
            isExist = os.path.exists(pathSave)
            if not isExist:
                os.makedirs(pathSave)

            # save Figures
            trfig.savefig(
                pathSave + str(trfig.get_label()) + "_position" + str(position) + ".png"
            )
            phfig.savefig(pathSave + str(phfig.get_label()) + ".png")
            sdfig.savefig(pathSave + str(sdfig.get_label()) + ".png")
            g2fig.savefig(pathSave + str(g2fig.get_label()) + ".png")

            # DATA ALLOCATION IN A DATA FRAME

            # We add the relevant parameters to a data frame
            efficiencyDF["# lcH"], efficiencyDF["# lcV"], efficiencyDF["# lc"] = (
                [sum(photonlcH)],
                sum(photonlcV),
                sum(photonlcH) + sum(photonlcV),
            )
            efficiencyDF["# kcH"], efficiencyDF["# kcV"], efficiencyDF["# kc"] = (
                sum(photonkcH),
                sum(photonkcV),
                sum(photonkcH) + sum(photonkcV),
            )
            efficiencyDF["eff lc"], efficiencyDF["std lc"] = photonlcEff, photonlcEfferr
            efficiencyDF["eff kc"], efficiencyDF["std kc"] = photonkcEff, photonkcEfferr
            # efficiencyDF['# (kc&lc)'] = coinc_n
            # efficiencyDF['eff (lc|kc)'], efficiencyDF['std (lc|kc)'] = eff_lc_upon_kc, eff_lc_upon_kcerr
            # efficiencyDF['eff (kc|lc)'], efficiencyDF['std (kc|lc)'] = eff_kc_upon_lc, eff_kc_upon_lcerr
            # efficiencyDF['eff (lc&kc)'], efficiencyDF['std (lc&kc)'] = coinc_n/fstot, np.sqrt(coinc_n)/fstot

            sdDF["sd triggers"] = [sdTrigtot]

            g2DF["# (lc>1) t=0"], g2DF["# (kc>1) t=0"] = (
                [coincidences_lc[10]],
                coincidences_kc[10],
            )
            g2DF["g2lc t=0"], g2DF["std g2lc t=0"] = g2_lc[10], errg2_lc[10]
            g2DF["g2kc t=0"], g2DF["std g2kc t=0"] = g2_kc[10], errg2_kc[10]

            # We save the data frame to an excelfile
            writer = pd.ExcelWriter(
                pathSave + filename + "position" + str(position) + "_results.xlsx"
            )
            efficiencyDF.to_excel(writer, sheet_name="efficiency")
            sdDF.to_excel(writer, sheet_name="SD")
            g2DF.to_excel(writer, sheet_name="g2")
            writer._save()

        # plt.show()
        if heralded:
            return {
                "F1, H": kclcHF1,
                "F1, V": kclcVF1,
                "F2, H": kclcHF2,
                "F2, V": kclcVF2,
            }
        else:
            return {"F1, H": lcHF1, "F1, V": lcVF1, "F2, H": lcHF2, "F2, V": lcVF2}

    def dataEv_MWRamseySequence(self, path: str, filename: str, filetype: str = "h5"):
        """
        INPUTS
        path: path to directory where to find the dataDic
        filename: name of the file to be analyzed
        filetype: h5 default

        OUTPUTS
        MWdurations: list of all the microwave pulse duration
        SDmean: list with all the mean state detection outcomes
        err_SDmean: list with errors for SDmean
        """
        # ------ Definition of parameters ------

        file_postSelected = (
            path + "\\goodAtomSelectorFiles\\" + filename + "_goodAtoms.pkl"
        )
        if not os.path.exists(file_postSelected):
            self.dataEv_postSelection(path, filename, filetype, kcCounts=2000, no=30)

        a_file = open(file_postSelected, "rb")
        atomDic = pickle.load(a_file)

        # Sequence timings
        coolingDur = (
            1.0550e-3 + 31.8 * 1e-6 - 1 * 1e-6
        )  # (s) duration of the atom cooling stage
        OPdur = 200e-6
        SDdur = 7.5e-6  # duration of the state detection stage
        endGap = 2e-6

        timeIntervalStart = 1086 * 1e-6  # (s) initial microwave duration
        NrPoints = 29  # microwave duration steps per ramp
        timeIntervalStep = 1e-6  # (s) microwave duration step

        # ------ Simple calculations of parameters ------
        timeIntervals = [
            timeIntervalStart + i * timeIntervalStep for i in range(0, NrPoints)
        ]
        SDstart = [coolingDur + OPdur + m for m in timeIntervals]
        SDend = [start + SDdur for start in SDstart]
        seqDur = SDend[-1] + endGap
        # Initilization of state detection arrays: each element corresponds to a trial, 0=dark, 1=bright
        SD_binary = dict((i, []) for i in range(NrPoints))
        SDthreshold = (
            1  # number of photons starting from which a bright state is witnessed
        )

        # We load the data in a dictionary - dataDic = {'label': [channel,[timeStamp]]}
        load = DataLoading()
        dataDic = load.data_goodAtoms(path, filename, atomDic)

        # dataDic = gf.data_goodAtoms(path, filename, atomDic)

        # We plot the photon detection time histogram
        syncRamp = self.syncFast2
        trigger = self.syncFast

        maxTrigDiff = seqDur + 100e-6  # should be bigger than triggers time difference
        binsize = 0.5 * 1e-6
        binNum = int(maxTrigDiff / binsize)
        detectors = [self.kcH, self.kcV, self.lcH, self.lcV]
        colors = ["violet", "violet", "tab:blue", "tab:blue"]
        fsdelay = {self.kcH: 0, self.kcV: 12e-9, self.lcH: 0, self.lcV: 0.0}

        gates = SDstart + SDend
        self.channels_histo(
            dataDic,
            detectors,
            gates,
            binNum,
            trigger,
            maxTrigDiff,
            fsdelay,
            filename,
            colors,
        )

        # ------ State detection ------

        # Loop over different scan triggers
        for sfs in tqdm(dataDic[syncRamp][1][:-1], file=sys.stdout):
            # start of each spectrum
            indexSF0 = np.searchsorted(dataDic[self.syncFast][1], sfs)
            for step in range(0, NrPoints):
                indexSF = indexSF0 + step
                if (
                    dataDic[self.syncFast][1][indexSF]
                    - dataDic[self.syncFast][1][indexSF - 1]
                    > SDend[-1] + 100e-6
                ):
                    print(
                        dataDic[self.syncFast][1][indexSF]
                        - dataDic[self.syncFast][1][indexSF - 1]
                    )
                    print(step)
                    print("BROKEN")
                    break

                sf = dataDic[self.syncFast][1][indexSF]

                left_h, right_h = np.searchsorted(
                    dataDic[self.kcH][1],
                    [
                        sf + fsdelay[self.kcH] + SDstart[step],
                        sf + fsdelay[self.kcH] + SDend[step],
                    ],
                )
                left_v, right_v = np.searchsorted(
                    dataDic[self.kcV][1],
                    [
                        sf + fsdelay[self.kcV] + SDstart[step],
                        sf + fsdelay[self.kcV] + SDend[step],
                    ],
                )

                if (right_h - left_h + right_v - left_v) >= SDthreshold:
                    SD_binary[step] = np.append(SD_binary[step], 1)
                else:
                    SD_binary[step] = np.append(SD_binary[step], 0)

        SDmean = np.zeros(NrPoints)
        err_SDmean = np.zeros(NrPoints)
        for k in range(0, NrPoints, 1):
            SDmean[k] = np.mean(SD_binary[k])
            err_SDmean[k] = np.std(SD_binary[k]) / np.sqrt(len(SD_binary[k]))

        print([sum(SD_binary[i]) for i in range(len(SD_binary))])
        print(SDmean)

        return timeIntervals, SDmean, err_SDmean

    def dataEv_MWspectrum(self, path: str, filename: str, filetype: str = "h5"):
        """
        # INPUTS
        # path: directory where to find the dataDic
        # filename: name of the file to be analyzed
        # filetype
        # MWdur: duration of the microwave pulse
        # MWfreqStepsPerRamp: number of microwave frequencies per ramp
        # MWfreqStep: frequency increment of the microwave
        # OUTPUTS
        # SDmean: list with all the mean state detection outcomes
        # err_SDmean: list with errors for SDmean
        """
        # ------ Definition of parameters ------

        file_postSelected = (
            path + "\\goodAtomSelectorFiles\\" + filename + "_goodAtoms.pkl"
        )
        if not os.path.exists(file_postSelected):
            self.dataEv_postSelection(path, filename, filetype, kcCounts=2500, no=100)

        a_file = open(file_postSelected, "rb")
        atomDic = pickle.load(a_file)
        # ------ Definition of parameter values ------
        MWfreqCenter = 0  # in kHz
        # MWfreqSpan = 300  # in kHz
        MWfreqSpan = 400  # in kHz
        MWfreqStart = MWfreqCenter - MWfreqSpan / 2  # (KHz) initial microwave duration
        MWfreqEnd = MWfreqCenter + MWfreqSpan / 2  # (KHz) initial microwave duration
        # MWfreqStepsPerRamp = 600  # number of microwave  steps per ramp
        MWfreqStepsPerRamp = 800  # number of microwave  steps per ramp
        MWfreqStep = MWfreqSpan / MWfreqStepsPerRamp  # (kHz) microwave step

        SDthreshold = (
            1  # number of photons starting from which a bright state is witnessed
        )

        # Sequence timings
        coolingDur = 400e-6
        OPdur = 100e-6
        gap1 = 1e-6
        gap2 = 3e-6
        endGap = 1e-6
        MWdur = 450e-6
        SDdur = 7.5e-6  # duration of the state detection stage
        seqDur = coolingDur + gap1 + OPdur + MWdur + gap2 + SDdur + endGap

        # ------ Simple calculations of parameters ------
        # Calculation of the microwave frequencies
        MWfrequencies = (
            np.linspace(MWfreqStart, MWfreqEnd, MWfreqStepsPerRamp) - MWfreqCenter
        )
        # Calculation of the state detection windows
        SDstart = seqDur - SDdur - endGap
        SDend = SDstart + SDdur

        print("SeqDur = %.2f us " % (seqDur * 1e6))
        print("SDstart = %.2f us " % (SDstart * 1e6))
        print("SDend = %.2f us " % (SDend * 1e6))
        print("MW step = %.3f kHz" % (MWfreqStep))

        # Initilization of state detection arrays: each element corresponds to a trial, 0=dark, 1=bright
        SD_binary = dict((i, []) for i in range(MWfreqStepsPerRamp))

        specDuration = seqDur * MWfreqStepsPerRamp  # s
        print("Duration 1 Spectrum = %.3f s" % (specDuration))
        atomDic = self.load.atomDicCleaner(atomDic, specDuration)
        print(atomDic)

        # We load the data in a dictionary - dataDic = {'label': [channel,[timeStamp]]}
        dataDic = self.load.data_goodAtoms(path, filename, atomDic)

        # We plot the photon detection time histogram
        syncRamp = self.syncFast2
        trigger = self.syncFast

        maxTrigDiff = seqDur + 100e-6  # should be bigger than triggers time difference
        binsize = 0.5 * 1e-6
        binNum = int(maxTrigDiff / binsize)
        detectors = [self.kcH, self.kcV, self.lcH, self.lcV]
        colors = ["violet", "violet", "tab:blue", "tab:blue"]
        fsdelay = {self.kcH: 0, self.kcV: 12e-9, self.lcH: 0, self.lcV: 0.0}

        gates = [SDstart, SDend]
        self.channels_histo(
            dataDic,
            detectors,
            gates,
            binNum,
            trigger,
            maxTrigDiff,
            fsdelay,
            filename,
            colors,
        )
        # plt.show()

        # ------ State detection ------

        # Loop over different scan triggers
        for sfs in tqdm(dataDic[syncRamp][1][:-2], file=sys.stdout):
            # start of each spectrum
            indexSF0 = np.searchsorted(dataDic[trigger][1], sfs)
            for step in range(0, MWfreqStepsPerRamp, 1):
                # iterate through all steps in a spectrum

                indexSF = indexSF0 + step * 1
                left_h = np.searchsorted(
                    dataDic[self.kcH][1],
                    dataDic[trigger][1][indexSF] + SDstart + fsdelay[self.kcH],
                )
                right_h = np.searchsorted(
                    dataDic[self.kcH][1],
                    dataDic[trigger][1][indexSF] + SDend + fsdelay[self.kcH],
                )
                left_v = np.searchsorted(
                    dataDic[self.kcV][1],
                    dataDic[trigger][1][indexSF] + SDstart + fsdelay[self.kcV],
                )
                right_v = np.searchsorted(
                    dataDic[self.kcV][1],
                    dataDic[trigger][1][indexSF] + SDend + fsdelay[self.kcV],
                )

                if right_h - left_h + right_v - left_v >= SDthreshold:
                    SD_binary[step] = np.append(SD_binary[step], 1)
                else:
                    SD_binary[step] = np.append(SD_binary[step], 0)

        SDmean = np.zeros(MWfreqStepsPerRamp)
        err_SDmean = np.zeros(MWfreqStepsPerRamp)
        for k in range(0, MWfreqStepsPerRamp, 1):
            SDmean[k] = np.mean(SD_binary[k])
            err_SDmean[k] = np.std(SD_binary[k]) / np.sqrt(len(SD_binary[k]))

        return MWfrequencies, SDmean, err_SDmean

    def dataEv_writeAtom(
        self, path: str, file_list: List[str], filetype: str = "h5", position: int = 0
    ):
        """
        INPUTS
        path: path to directory where to find the dataDic
        filename: name of the file to be analyzed
        filetype: h5 default

        OUTPUT:
        Whatever
        """
        pathSave = path + "\\atomStateFiles\\"
        efficiencyDF = pd.DataFrame()
        sdDF = pd.DataFrame()
        g2DF = pd.DataFrame()

        print("\nAnalysing: ")
        for filename in file_list:
            print(filename + ", ")
        print("trigger ", position + 1)

        for filename in file_list:
            file_postSelected = (
                path + "\\goodAtomSelectorFiles\\" + filename + "_goodAtoms.pkl"
            )
            if not os.path.exists(file_postSelected):
                self.dataEv_postSelection(path, filename, filetype)

            a_file = open(file_postSelected, "rb")
            atomDic = pickle.load(a_file)
            print(atomDic)

            # ------ We get the data ------
            dataDic = self.load.data_goodAtoms(path, filename, atomDic)

            cooldur = 1.0034e-3
            OPdur = 150e-6
            pulseDur = 100e-9
            SDDur = 7.5e-6

            opgate = [cooldur, cooldur + OPdur]
            writegate = [opgate[1], opgate[1] + pulseDur]
            SDgate = [1.2384e-3, 1.2384e-3 + SDDur]
            coincgate = [1.15430e-3, 1.1548e-3]

            # gates = self.flatten(writegate) #+ self.flatten(SDgate) + self.flatten(coincgate)
            gates = writegate + SDgate + coincgate

            seqDur = gates[-1] + 100e-6
            seqDur = 1.3e-3 + 100e-6
            print(seqDur)

            trigger = self.syncFast
            maxTrigDiff = seqDur  # should be bigger than triggers time difference
            binsize = 0.01 * 1e-6
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
            chfig = self.channels_histo(
                dataDic,
                detectors,
                gates,
                binNum,
                trigger,
                maxTrigDiff,
                fsdelay,
                filename,
                colors,
            )
            # plt.show(block = True)

            fstot = 0
            sdTrigtot = 0
            tlcH, tlcV, tkcH, tkcV = [], [], [], []
            photonlcH, photonlcV, photonkcH, photonkcV = [], [], [], []
            photonSD = [[], []]  # n of photons in SD gate for each trial
            kcHF1, kcHF2, kcVF1, kcVF2 = 0, 0, 0, 0
            SDthreshold = 1

            print("""
                ============
                Atoms data
                ============
    
                --- Looping over Fast Sequence triggers ---     
                """)

            # Loop over all the fast sequence triggers
            # for fsr in tqdm(dataDic[trig][1][:-5], file=sys.stdout):
            #     fs = dataDic[self.syncFast][1][
            #         np.searchsorted(dataDic[self.syncFast][1], fsr) + position]  # + position*1e-6

            trig = self.syncFast2
            # trig = self.syncFast

            time.sleep(0.1)
            # Loop over all the fast sequence triggers
            for fsr in tqdm(dataDic[trig][1][:-5], file=sys.stdout):
                fs = dataDic[self.syncFast][1][
                    np.searchsorted(dataDic[self.syncFast][1], fsr) + position
                ]  # + position*1e-6
                # kcH herald
                sf = fs + fsdelay[self.kcH]
                left, right = np.searchsorted(
                    dataDic[self.kcH][1], [sf + writegate[0], sf + writegate[1]]
                )
                # (tkcH.extend(dataDic[self.kcH][1][left:right] - sf) if left - right != 0 else None)
                photonkcH.append(right - left)
                # kcV herald
                sf = fs + fsdelay[self.kcV]
                left, right = np.searchsorted(
                    dataDic[self.kcV][1], [sf + writegate[0], sf + writegate[1]]
                )
                # (tkcV.extend(dataDic[self.kcV][1][left:right] - sf) if left - right != 0 else None)
                photonkcV.append(right - left)
                # SD
                left_h, right_h = np.searchsorted(
                    dataDic[self.kcH][1], [sf + SDgate[0], sf + SDgate[1]]
                )
                left_v, right_v = np.searchsorted(
                    dataDic[self.kcV][1], [sf + SDgate[0], sf + SDgate[1]]
                )
                left_cond, right_cond = np.searchsorted(
                    dataDic[self.sdTrig][1], [sf + coincgate[0], sf + coincgate[1]]
                )

                if right_cond - left_cond >= 1:
                    sdTrigtot += 1

                    photonSD[0].append((right_v - left_v) + (right_h - left_h))
                    (
                        photonSD[1].extend(dataDic[self.kcH][1][left_h:right_h] - sf)
                        if left_h - right_h != 0
                        else None
                    )
                    (
                        photonSD[1].extend(dataDic[self.kcV][1][left_v:right_v] - sf)
                        if left_v - right_v != 0
                        else None
                    )

                    # heralded coincidences
                    kcHF1 += photonkcH[-1] & (photonSD[0][-1] == 0)
                    kcHF2 += photonkcH[-1] & (photonSD[0][-1] >= SDthreshold)
                    kcVF1 += photonkcV[-1] & (photonSD[0][-1] == 0)
                    kcVF2 += photonkcV[-1] & (photonSD[0][-1] >= SDthreshold)

                fstot += 1

        effkcH, effkcH_err = sum(photonkcH) / fstot, np.sqrt(sum(photonkcH)) / fstot
        effkcV, effkcV_err = sum(photonkcV) / fstot, np.sqrt(sum(photonkcV)) / fstot
        photonkcEff, photonkcEfferr = (
            (sum(photonkcV) + sum(photonkcH)) / fstot,
            np.sqrt(sum(photonkcV) + sum(photonkcH)) / fstot,
        )

        print("Efficiency kcH = %.4f(%d)" % (effkcH, effkcH_err * 1e4))
        print("Efficiency kcV = %.4f(%d)" % (effkcV, effkcV_err * 1e4))
        print(
            "Efficiency herald total = %.4f(%d)" % (photonkcEff, photonkcEfferr * 1e4)
        )

        pF2, pF2err = (
            (kcHF2 + kcVF2) / (kcHF1 + kcHF2 + kcVF1 + kcVF2),
            np.sqrt(kcHF2 + kcVF2) / (kcHF1 + kcHF2 + kcVF1 + kcVF2),
        )
        pF1, pF1err = (
            (kcHF1 + kcVF1) / (kcHF1 + kcHF2 + kcVF1 + kcVF2),
            np.sqrt(kcHF1 + kcVF1) / (kcHF1 + kcHF2 + kcVF1 + kcVF2),
        )

        print("Probability F1 = %.3f(%d)" % (pF1, pF1err * 1e3))
        print("Probability F2 = %.3f(%d)" % (pF2, pF2err * 1e3))

        # --- g2 analysis ---
        g2zero, g2zeroerr, g2fig = self.g2(photonkcH, photonkcV, seqDur, filename)

        if self.dataSave:
            # Check whether the specified path exists or not
            isExist = os.path.exists(pathSave)
            if not isExist:
                os.makedirs(pathSave)

            # save Figures
            chfig.savefig(pathSave + str(chfig.get_label()) + ".png")
            g2fig.savefig(pathSave + str(g2fig.get_label()) + ".png")

            # DATA ALLOCATION IN A DATA FRAME

            # We add the relevant parameters to a data frame
            efficiencyDF["# kcH"], efficiencyDF["# kcV"], efficiencyDF["# kc"] = (
                [sum(photonkcH)],
                sum(photonkcV),
                sum(photonkcH) + sum(photonkcV),
            )
            efficiencyDF["eff kc"], efficiencyDF["std kc"] = photonkcEff, photonkcEfferr
            efficiencyDF["prob kc&F2"], efficiencyDF["std kc&F2"] = pF2, pF2err
            efficiencyDF["prob kc&F1"], efficiencyDF["std kc&F1"] = pF1, pF1err

            sdDF["sd triggers"] = [sdTrigtot]

            # g2DF['# (kc>1) t=0'] = [coincidences_10]
            g2DF["g2 t=0"], g2DF["std g2 t=0"] = [g2zero], g2zeroerr
            # We save the data frame to an excelfile
            writer = pd.ExcelWriter(pathSave + filename + "_results.xlsx")
            efficiencyDF.to_excel(writer, sheet_name="efficiency")
            sdDF.to_excel(writer, sheet_name="SD")
            g2DF.to_excel(writer, sheet_name="g2")
            writer._save()

        return {"F1": kcHF1 + kcVF1, "F2": kcHF2 + kcVF2}

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
        # print(diffFS)
        maxFSdur = np.amax(
            diffFS[diffFS < maxTrigDiff]
        )  # time difference between atoms are excluded
        # maxFSdur = 45e-3
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

    def g2(self, counts1, counts2, seqDur, filename: str = " "):
        # ---- g2 autocorrelation analysis ----
        trialShifts = np.arange(-10, 10)
        coincidences = np.zeros(len(trialShifts))

        for i, shift in enumerate(trialShifts):
            coincidences[i] = np.inner(counts1, np.roll(counts2, shift))

        g2 = coincidences / np.average(np.delete(coincidences, 10))
        g2err = (
            np.clip(coincidences, 1, 1e15)
            / np.average(np.delete(coincidences, 10))
            * np.sqrt(
                1 / np.clip(coincidences, 1, 1e15)
                + 1 / np.sum(np.delete(coincidences, 10))
            )
        )

        # print("""
        #        ============
        #        g2 Analysis
        #        ============
        #        """)
        #
        # print("Photons in kcPi = %d\nPhotons in kcV = %d" % (sum(photonkcH), sum(photonkcV)))
        # print("coincidences(0) kC = ", coincidences_kc[10])
        text = "g2(0) = %.2f(%d)" % (g2[10], g2err[10] * 1e2)
        print(text)

        spttl = filename + " - g2"
        fig = plt.figure(spttl, figsize=[10, 7])
        ax = fig.add_subplot(1, 1, 1)

        fig.suptitle(spttl + "\n" + text)
        ax.bar(trialShifts, g2, yerr=g2err, color=self.colour["orangeDark"])
        ax.set_xlabel(r"$\tau$ (ms)")

        plt.tight_layout()

        return g2[10], g2err[10], fig

    # --- Fitting procedures ---
    def dataFit_exponential(self, x, y, yerr):

        fit = fittingFunctions()
        # Exponential fit
        model = lmfit.Model(fit.exponential)  # model used for fitting
        params = model.make_params()

        # Fitting parameter limits definition
        def pset(name, value, min=None, max=None, vary=True):
            params[name].set(value=value, min=min, max=max, vary=vary)

        # Parameters for the exponential fit pset(label,start,lowerLimit,upperLimit, fittingBoolean)
        pset("A0", 0.8, 0.2, 1, True)
        pset("tau", 0.2, 0.1, 1, True)
        pset("offset", 0.5, 0, 0.6, True)

        results = model.fit(y, params, x=x, weights=1 / yerr**2, options={"maxfev": 1})
        print(results.fit_report())
        return (model, results)

    def dataFit_rabiOsc(self, x, y):
        # INPUTS
        # MWdurations: list of the microwave duration ramp values
        # P0: probability that the atom is in state 1 (F=2)
        # OUTPUTS
        # modRabiOsc: fit model
        # res_rabiOsc: fit result
        fit = fittingFunctions()
        model = lmfit.Model(fit.rabiRotation)  # model used for fitting
        params = model.make_params()

        # Fitting parameter limits definition
        def pset(name, value, min=None, max=None, vary=True):
            params[name].set(value=value, min=min, max=max, vary=vary)

        # Parameters for the rabi Oscillation fit (label,start,lowerLimit,upperLimit, fittingBoolean)
        pset("freq", 7, 1, 200, True)
        pset("decRate", 1, 0, 500, True)
        pset("A", 0.92, 0.001, 1, True)
        pset("off", 0.00, 0.0, 0.1, False)

        print(params)
        results = model.fit(y, params, time=x)
        print(results.fit_report())

        return (model, results)

    def dataFit_ramseyOsc(self, x, y):
        # INPUTS
        # MWdurations: list of the microwave duration ramp values
        # P0: probability that the atom is in state 1 (F=2)
        # OUTPUTS
        # modRabiOsc: fit model
        # res_rabiOsc: fit result
        fit = fittingFunctions()
        model = lmfit.Model(fit.ramseyRotation)  # model used for fitting
        params = model.make_params()

        # Fitting parameter limits definition
        def pset(name, value, min=None, max=None, vary=True):
            params[name].set(value=value, min=min, max=max, vary=vary)

        # Parameters for the rabi Oscillation fit (label,start,lowerLimit,upperLimit, fittingBoolean)
        pset("freq", 25.0, 1, 200, True)
        pset("decRate", 1.3, 0, 500, True)
        pset("A", 0.95, 0.001, 2, True)
        pset("off", 0.0, -0.1, 0.1, False)

        print(params)
        results = model.fit(y, params, time=x)
        print(results.fit_report())

        return (model, results)

    def dataFit_lorenzian(self, x, y, yerr):

        # Lorenzian fit
        fit = fittingFunctions()
        model = lmfit.Model(fit.lorenzian)  # model used for fitting
        params = model.make_params()

        # Fitting parameter limits definition
        def pset(name, value, min=None, max=None, vary=True):
            params[name].set(value=value, min=min, max=max, vary=vary)

        # Parameters for the cavity locked fit pset(label,start,lowerLimit,upperLimit, fittingBoolean)
        pset("amp1", 0.35, 0.05, 3, True)  # arbitrary amplitude
        pset("k1", 1.5, 0.1, 30, True)  # total cavity decay rate
        pset("offset", 0.01, 0, 0.5, True)  # total cavity decay rate

        l = min(x)
        r = max(x)
        c = (l + r) / 2

        pset("foff", c, l, r, True)

        # results = model.fit(y,params,freq=x, weights=1/yerr**2,options={'maxfev': 1})
        results = model.fit(y, params, freq=x, options={"maxfev": 1})
        print(results.fit_report())
        return (model, results)

    # --- Useful functions ---
    def flatten(self, xss):
        return [x for xs in xss for x in xs]

    def find_local_maxima_in_range(self, data, window_size, threshold):
        maxima_indices = []
        maxima_values = []

        # Process the data in non-overlapping windows
        for i in range(0, len(data), window_size):
            print(i)
            # Extract the current window from the data
            window = data[i : i + window_size]

            # Find the maximum value in the current window
            max_value = np.max(window)

            # Check if the maximum value is above the threshold
            if max_value > threshold:
                # Get the index of the maximum value within the window
                max_index_in_window = np.argmax(window)

                # Convert the index to the global index in the original data
                global_index = i + max_index_in_window

                # Store the global index and the maximum value
                maxima_indices.append(global_index)
                maxima_values.append(max_value)

        return np.array(maxima_indices), np.array(maxima_values)


class stateTomography:
    def __init__(self):

        # ------ Definition of parameters and useful matrices------
        self.pauli = [qeye(2), sigmax(), sigmay(), sigmaz()]

    def stokes_1q(self, counts: dict, prntP: bool):

        # from counts, we compute the probabilities with errors
        p = {}
        for base in counts.keys():
            p[base] = {"F1": {"p": 0, "err": 0}, "F2": {"p": 0, "err": 0}}

        for base in counts.keys():
            for state in p[base].keys():
                p[base][state]["p"] = counts[base][state] / sum(
                    [n for n in counts[base].values()]
                )
                p[base][state]["err"] = np.sqrt(counts[base][state]) / sum(
                    [n for n in counts[base].values()]
                )

        # We compute the  Stokes parameters from the data we have collected
        S = [0 for x in range(4)]  # Stokes parameters for 1 qubit state
        S_err = [0 for x in range(4)]  # Stokes parameters errors for 1 qubit state

        S[0] = p["Z"]["F2"]["p"] + p["Z"]["F1"]["p"]

        S[1] = p["X"]["F2"]["p"] - p["X"]["F1"]["p"]

        S[2] = p["Y"]["F2"]["p"] - p["Y"]["F1"]["p"]

        S[3] = p["Z"]["F2"]["p"] - p["Z"]["F1"]["p"]

        S_err[0] = np.sqrt(p["Z"]["F1"]["err"] ** 2 + p["Z"]["F2"]["err"] ** 2)
        S_err[1] = np.sqrt(p["X"]["F1"]["err"] ** 2 + p["X"]["F2"]["err"] ** 2)
        S_err[2] = np.sqrt(p["Y"]["F1"]["err"] ** 2 + p["Y"]["F2"]["err"] ** 2)
        S_err[3] = np.sqrt(p["Z"]["F1"]["err"] ** 2 + p["Z"]["F2"]["err"] ** 2)

        if prntP:
            for x in p:
                print(x)
                for y in p[x]:
                    print(y, ":", p[x][y])
        return S, S_err

    def stokes_2q(self, coincidences: dict, prntP: bool):

        # from coincidences, we compute the probabilities with errors
        p = {}
        for base in coincidences.keys():
            p[base] = {
                "F1, H": {"p": 0, "err": 0},
                "F1, V": {"p": 0, "err": 0},
                "F2, H": {"p": 0, "err": 0},
                "F2, V": {"p": 0, "err": 0},
            }
        for base in coincidences.keys():
            for count in p[base].keys():
                p[base][count]["p"] = coincidences[base][count] / sum(
                    [n for n in coincidences[base].values()]
                )
                p[base][count]["err"] = np.sqrt(coincidences[base][count]) / sum(
                    [n for n in coincidences[base].values()]
                )

        # We compute the  Stokes parameters from the data we have collected
        S = [[0 for x in range(4)] for y in range(4)]  # Stokes parameters
        S_err = [[0 for x in range(4)] for y in range(4)]  # Stokes parameters errors
        S[0][0] = 1

        S[1][1] = (
            p["X_DA"]["F1, H"]["p"]
            - p["X_DA"]["F1, V"]["p"]
            - p["X_DA"]["F2, H"]["p"]
            + p["X_DA"]["F2, V"]["p"]
        )
        S[1][2] = (
            p["X_HV"]["F1, H"]["p"]
            - p["X_HV"]["F1, V"]["p"]
            - p["X_HV"]["F2, H"]["p"]
            + p["X_HV"]["F2, V"]["p"]
        )
        S[1][3] = (
            p["X_RL"]["F1, H"]["p"]
            - p["X_RL"]["F1, V"]["p"]
            - p["X_RL"]["F2, H"]["p"]
            + p["X_RL"]["F2, V"]["p"]
        )

        S[2][1] = (
            p["Y_DA"]["F1, H"]["p"]
            - p["Y_DA"]["F1, V"]["p"]
            - p["Y_DA"]["F2, H"]["p"]
            + p["Y_DA"]["F2, V"]["p"]
        )
        S[2][2] = (
            p["Y_HV"]["F1, H"]["p"]
            - p["Y_HV"]["F1, V"]["p"]
            - p["Y_HV"]["F2, H"]["p"]
            + p["Y_HV"]["F2, V"]["p"]
        )
        S[2][3] = (
            p["Y_RL"]["F1, H"]["p"]
            - p["Y_RL"]["F1, V"]["p"]
            - p["Y_RL"]["F2, H"]["p"]
            + p["Y_RL"]["F2, V"]["p"]
        )

        S[3][1] = (
            p["Z_DA"]["F1, H"]["p"]
            - p["Z_DA"]["F1, V"]["p"]
            - p["Z_DA"]["F2, H"]["p"]
            + p["Z_DA"]["F2, V"]["p"]
        )
        S[3][3] = (
            p["Z_RL"]["F1, H"]["p"]
            - p["Z_RL"]["F1, V"]["p"]
            - p["Z_RL"]["F2, H"]["p"]
            + p["Z_RL"]["F2, V"]["p"]
        )
        S[3][2] = (
            p["Z_HV"]["F1, H"]["p"]
            - p["Z_HV"]["F1, V"]["p"]
            - p["Z_HV"]["F2, H"]["p"]
            + p["Z_HV"]["F2, V"]["p"]
        )

        S[0][1] = (
            p["X_DA"]["F1, H"]["p"]
            - p["X_DA"]["F1, V"]["p"]
            + p["X_DA"]["F2, H"]["p"]
            - p["X_DA"]["F2, V"]["p"]
        )
        S[0][2] = (
            p["Y_HV"]["F1, H"]["p"]
            - p["Y_HV"]["F2, H"]["p"]
            + p["Y_HV"]["F1, V"]["p"]
            - p["Y_HV"]["F2, V"]["p"]
        )
        S[0][3] = (
            p["Z_RL"]["F1, H"]["p"]
            - p["Z_RL"]["F1, V"]["p"]
            + p["Z_RL"]["F2, H"]["p"]
            - p["Z_RL"]["F2, V"]["p"]
        )
        S[1][0] = (
            p["X_DA"]["F1, H"]["p"]
            + p["X_DA"]["F1, V"]["p"]
            - p["X_DA"]["F2, H"]["p"]
            - p["X_DA"]["F2, V"]["p"]
        )
        S[2][0] = (
            p["Y_HV"]["F1, H"]["p"]
            + p["Y_HV"]["F2, H"]["p"]
            - p["Y_HV"]["F1, V"]["p"]
            - p["Y_HV"]["F2, V"]["p"]
        )
        S[3][0] = (
            p["Z_RL"]["F1, H"]["p"]
            + p["Z_RL"]["F1, V"]["p"]
            - p["Z_RL"]["F2, H"]["p"]
            - p["Z_RL"]["F2, V"]["p"]
        )

        S_err[1][1] = np.sqrt(
            p["X_DA"]["F1, H"]["err"] ** 2
            + p["X_DA"]["F1, V"]["err"] ** 2
            + p["X_DA"]["F2, H"]["err"] ** 2
            + p["X_DA"]["F2, V"]["err"] ** 2
        )
        S_err[1][2] = np.sqrt(
            p["X_HV"]["F1, H"]["err"] ** 2
            + p["X_HV"]["F1, V"]["err"] ** 2
            + p["X_HV"]["F2, H"]["err"] ** 2
            + p["X_HV"]["F2, V"]["err"] ** 2
        )
        S_err[1][3] = np.sqrt(
            p["X_RL"]["F1, H"]["err"] ** 2
            + p["X_RL"]["F1, V"]["err"] ** 2
            + p["X_RL"]["F2, H"]["err"] ** 2
            + p["X_RL"]["F2, V"]["err"] ** 2
        )

        S_err[2][1] = np.sqrt(
            p["Y_DA"]["F1, H"]["err"] ** 2
            + p["Y_DA"]["F1, V"]["err"] ** 2
            + p["Y_DA"]["F2, H"]["err"] ** 2
            + p["Y_DA"]["F2, V"]["err"] ** 2
        )
        S_err[2][2] = np.sqrt(
            p["Y_HV"]["F1, H"]["err"] ** 2
            + p["Y_HV"]["F1, V"]["err"] ** 2
            + p["Y_HV"]["F2, H"]["err"] ** 2
            + p["Y_HV"]["F2, V"]["err"] ** 2
        )
        S_err[2][3] = np.sqrt(
            p["Y_RL"]["F1, H"]["err"] ** 2
            + p["Y_RL"]["F1, V"]["err"] ** 2
            + p["Y_RL"]["F2, H"]["err"] ** 2
            + p["Y_RL"]["F2, V"]["err"] ** 2
        )

        S_err[3][1] = np.sqrt(
            p["Z_DA"]["F1, H"]["err"] ** 2
            + p["Z_DA"]["F1, V"]["err"] ** 2
            + p["Z_DA"]["F2, H"]["err"] ** 2
            + p["Z_DA"]["F2, V"]["err"] ** 2
        )
        S_err[3][3] = np.sqrt(
            p["Z_RL"]["F1, H"]["err"] ** 2
            + p["Z_RL"]["F1, V"]["err"] ** 2
            + p["Z_RL"]["F2, H"]["err"] ** 2
            + p["Z_RL"]["F2, V"]["err"] ** 2
        )
        S_err[3][2] = np.sqrt(
            p["Z_HV"]["F1, H"]["err"] ** 2
            + p["Z_HV"]["F1, V"]["err"] ** 2
            + p["Z_HV"]["F2, H"]["err"] ** 2
            + p["Z_HV"]["F2, V"]["err"] ** 2
        )

        S_err[0][1] = np.sqrt(
            p["X_DA"]["F1, H"]["err"] ** 2
            + p["X_DA"]["F1, V"]["err"] ** 2
            + p["X_DA"]["F2, H"]["err"] ** 2
            + p["X_DA"]["F2, V"]["err"] ** 2
        )
        S_err[0][2] = np.sqrt(
            p["Y_HV"]["F1, H"]["err"] ** 2
            + p["Y_HV"]["F2, H"]["err"] ** 2
            + p["Y_HV"]["F1, V"]["err"] ** 2
            + p["Y_HV"]["F2, V"]["err"] ** 2
        )
        S_err[0][3] = np.sqrt(
            p["Z_RL"]["F1, H"]["err"] ** 2
            + p["Z_RL"]["F1, V"]["err"] ** 2
            + p["Z_RL"]["F2, H"]["err"] ** 2
            + p["Z_RL"]["F2, V"]["err"] ** 2
        )
        S_err[1][0] = np.sqrt(
            p["X_DA"]["F1, H"]["err"] ** 2
            + p["X_DA"]["F1, V"]["err"] ** 2
            + p["X_DA"]["F2, H"]["err"] ** 2
            + p["X_DA"]["F2, V"]["err"] ** 2
        )
        S_err[2][0] = np.sqrt(
            p["Y_HV"]["F1, H"]["err"] ** 2
            + p["Y_HV"]["F2, H"]["err"] ** 2
            + p["Y_HV"]["F1, V"]["err"] ** 2
            + p["Y_HV"]["F2, V"]["err"] ** 2
        )
        S_err[3][0] = np.sqrt(
            p["Z_RL"]["F1, H"]["err"] ** 2
            + p["Z_RL"]["F1, V"]["err"] ** 2
            + p["Z_RL"]["F2, H"]["err"] ** 2
            + p["Z_RL"]["F2, V"]["err"] ** 2
        )

        if prntP:
            for x in p:
                print(x)
                for y in p[x]:
                    print(y, ":", p[x][y])

        return S, S_err

    def stokesBell(self):
        S = [[0 for x in range(4)] for y in range(4)]  # Stokes parameters
        S[0][0] = 1

        S[1][1] = 1
        S[1][2] = 0
        S[1][3] = 0

        S[2][1] = 0
        S[2][2] = -1
        S[2][3] = 0

        S[3][1] = 0
        S[3][2] = 0
        S[3][3] = 1

        S[0][1] = 0
        S[0][2] = 0
        S[0][3] = 0
        S[1][0] = 0
        S[2][0] = 0
        S[3][0] = 0

        return S

    def densitymatrix_1q(self, counts: dict):
        """return density matrix give a coincidences dictionary"""
        print("""
        =====================================      
        Computation - measured density matrix
        =====================================
        """)
        self.s, self.s_err = self.stokes_1q(counts, False)
        density = 0 * self.pauli[0]
        density_err = 0 * self.pauli[0]
        for i, sigma in enumerate(self.pauli):
            density = density + 1 / 2 * self.s[i] * sigma
            density_err = density_err + 1 / 4 * self.s_err[i] ** 2 * sigma
        density_err = Qobj.sqrtm(density_err)

        np.set_printoptions(precision=3)
        print("dm from stokes: ", density)
        print("dm_err from stokes: ", density_err)
        print("Eigenvalues of my density matrix = ", density.eigenenergies())

        return density, density_err

    def densitymatrix_2q(self, coincidences: dict):
        """return density matrix give a coincidences dictionary"""
        print("""
        =====================================      
        Computation - measured density matrix
        =====================================
        """)
        self.s, self.s_err = self.stokes_2q(coincidences, False)
        density = 0 * tensor(self.pauli[0], self.pauli[0])
        density_err = 0 * tensor(self.pauli[0], self.pauli[0])
        for i, sigma1 in enumerate(self.pauli):
            for j, sigma2 in enumerate(self.pauli):
                density = density + 1 / 4 * self.s[i][j] * tensor(sigma1, sigma2)
                density_err = density_err + 1 / 16 * self.s_err[i][j] ** 2 * tensor(
                    sigma1, sigma2
                )
        density_err = Qobj.sqrtm(density_err)

        np.set_printoptions(precision=3)
        print("dm from stokes: ", density)
        print("dm_err from stokes: ", density_err)
        print("Eigenvalues of my density matrix = ", density.eigenenergies())

        return density, density_err

    # We compute the angles we have to rotate our density matrix to obtain the maximum fidelity.
    def R(self, theta, phi):
        return gates.rx(theta) * gates.ry(phi)

    def fidelity_1q(self, state, dm, error=False):

        f = np.abs(state.dag() * dm * state)
        f_err = 0

        if error:
            for i, sigma in enumerate(self.pauli):
                f_err = (
                    f_err
                    + 1
                    / 4
                    * self.s_err[i] ** 2
                    * np.abs(state.dag() * sigma * state) ** 2
                )
        f_err = np.sqrt(f_err)

        return f, f_err

    def fidelity_2q(self, state, dm, error=False):

        f = np.abs(state.dag() * dm * state)
        f_err = 0

        if error:
            for i, sigma1 in enumerate(self.pauli):
                for j, sigma2 in enumerate(self.pauli):
                    f_err = (
                        f_err
                        + 1
                        / 16
                        * self.s_err[i][j] ** 2
                        * np.abs(state.dag() * tensor(sigma1, sigma2) * state) ** 2
                    )
            f_err = np.sqrt(f_err)

        return f, f_err

    def bestfidelity_1q(self, state, dm, error=False):

        def minus_fidelity(x):
            dm_r = 0 * dm  # measured density matrix rotated

            for i, sigma in enumerate(self.pauli):
                dm_r = (
                    dm_r
                    + 1
                    / 2
                    * self.s[i]
                    * self.R(x[0], x[1])
                    * sigma
                    * self.R(x[0], x[1]).dag()
                )

            f_r = self.fidelity_1q(state, dm_r)[0]

            return -f_r

        f_best, f_err = 0, 0
        dm_best = 0 * self.pauli[0]  # measured density matrix rotated, best

        x0 = np.array([np.radians(-47), np.radians(-70)])
        # x0 = np.array([np.radians(0), np.radians(0)])
        res = minimize(
            minus_fidelity,
            x0,
            method="nelder-mead",
            options={"xatol": 1e-8, "disp": True},
        )
        print("Rotation angles Atom (x,y) = ", np.degrees(res.x))
        f_best = -minus_fidelity(res.x)

        if error:
            for i, sigma in enumerate(self.pauli):
                f_err = (
                    f_err
                    + 1
                    / 4
                    * self.s_err[i] ** 2
                    * np.abs(state.dag() * sigma * state) ** 2
                )
            f_err = np.sqrt(f_err)

        for i, sigma in enumerate(self.pauli):
            dm_best = (
                dm_best
                + 1
                / 2
                * self.s[i]
                * self.R(res.x[0], res.x[1])
                * sigma
                * self.R(res.x[0], res.x[1]).dag()
            )

        return f_best, f_err, dm_best

    def bestfidelity_2q(self, state, dm, error=False):

        def minus_fidelity(x):
            dm_r = 0 * dm  # measured density matrix rotated

            for i, sigma1 in enumerate(self.pauli):
                for j, sigma2 in enumerate(self.pauli):
                    dm_r = dm_r + 1 / 4 * self.s[i][j] * tensor(
                        self.R(x[0], x[1]) * sigma1 * self.R(x[0], x[1]).dag(),
                        self.R(x[2], x[3]) * sigma2 * self.R(x[2], x[3]).dag(),
                    )

            f_r = self.fidelity_2q(state, dm_r)[0]

            return -f_r

        f_best, f_err = 0, 0
        dm_best = 0 * tensor(
            self.pauli[0], self.pauli[0]
        )  # measured density matrix rotated, best

        x0 = np.array([0, 0, 0, 0])
        res = minimize(
            minus_fidelity,
            x0,
            method="nelder-mead",
            options={"xatol": 1e-8, "disp": True},
        )
        print("Rotation angles Atoms (x,y) and Photon (x,y) = ", np.degrees(res.x))
        f_best = -minus_fidelity(res.x)

        if error:
            for i, sigma1 in enumerate(self.pauli):
                for j, sigma2 in enumerate(self.pauli):
                    f_err = (
                        f_err
                        + 1
                        / 16
                        * self.s_err[i][j] ** 2
                        * np.abs(
                            state.dag()
                            * tensor(
                                self.R(res.x[0], res.x[1])
                                * sigma1
                                * self.R(res.x[0], res.x[1]).dag(),
                                self.R(res.x[2], res.x[3])
                                * sigma2
                                * self.R(res.x[2], res.x[3]).dag(),
                            )
                            * state
                        )
                        ** 2
                    )

            f_err = np.sqrt(f_err)

        for i, sigma1 in enumerate(self.pauli):
            for j, sigma2 in enumerate(self.pauli):
                dm_best = dm_best + 1 / 4 * self.s[i][j] * tensor(
                    self.R(res.x[0], res.x[1])
                    * sigma1
                    * self.R(res.x[0], res.x[1]).dag(),
                    self.R(res.x[2], res.x[3])
                    * sigma2
                    * self.R(res.x[2], res.x[3]).dag(),
                )

        return f_best, f_err, dm_best


if __name__ == "__main__":
    """
    This script analyzes Rabi Flopping curves.
    """

    # ------ Setting the print options for numpy ------
    np.set_printoptions(legacy="1.25")

    # ------ Setting flags------
    analyze: bool = True
    # ------ Definition of data sources and destinations------
    path: str = os.path.abspath(os.path.dirname(__file__))
    filetype: str = ".h5"

    # ----- Printing flags -----
    if analyze:
        print("analizing from h5 files")
    else:
        print("Not Analyzing")

    # ------ Begin Analysis ------
    analysis = AtomAnalysis()

    fileList = ["16_02_26_MWRamsey_0_0_1"]

    for filename in fileList:
        print(filename)

    if analyze:
        #     for filename in fileList:
        (timeIntervals, SDmean, err_SDmean) = analysis.dataEv_MWRamseySequence(
            path, filename, filetype
        )

        # ------ Data is added to a data frame ------
        df = pd.DataFrame()
        df["timeInterval"] = timeIntervals
        df["SDmean"] = SDmean
        df["err_SDmean"] = err_SDmean

        # ------ Data fitting ----
        model, fit_results = analysis.dataFit_ramseyOsc(
            df["timeInterval"], df["SDmean"]
        )

        # --- Data Plotting ---
        timeFitPlot = np.linspace(0, df["timeInterval"][:-1] + 500e-6, 1000)
        plt.rcParams.update({"font.size": 15})

        fig1 = plt.figure(str(filename), figsize=(8, 4), constrained_layout=True)
        report = fit_results.fit_report(show_correl=False)
        pipulse = 1 / fit_results.best_values["freq"] * 1 / 4  # us
        index = report.find("[[Variables]]")
        print(index)
        print(report[index + len("[[Variables]]") :])
        print()
        fig1.suptitle(
            filename
            + report[index + len("[[Variables]]") :]
            + "\n"
            + "Pi pulse: %.1f us" % (pipulse * 1e3)
        )
        ax1 = fig1.add_subplot(1, 1, 1)
        ax1.plot(
            timeFitPlot * 1e6,
            model.eval(time=timeFitPlot, params=fit_results.params),
            color=analysis.colour["blueLight"],
            zorder=-1,
        )
        ax1.errorbar(
            df["timeInterval"] * 1e6,
            df["SDmean"],
            yerr=df["err_SDmean"],
            marker="o",
            linestyle="",
            color=analysis.colour["blueDark"],
            capsize=3,
        )
        ax1.set_xlabel("Microwave pulse duration ($\mathrm{\mu}$s)")
        ax1.set_ylabel(r"$\mathrm{P}\left(\left|F=2\right\rangle\right)$")
        ax1.tick_params(labelsize=20, direction="inout")
        ax1.grid()
        ax1.set_ylim([0, 1])
        ax1.set_xlim([0, 1200])

        fig1.savefig(path + "\\" + fig1.get_label() + ".pdf")
        fig1.savefig(path + "\\" + fig1.get_label() + ".png")

        plt.show(block=True)
