import sys
import h5py
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pprint import pprint


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

        filedata = h5py.File(path + fileName + ".h5", mode="r")

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
            dataD = self.data_loading(path, fileName, atom)

            dataD["ch1"][1] = np.array([x for x in dataD["ch1"][1]])

            # We take the time-stamps from the valid atom time period
            # index of the first syncFastScan that we don't consider
            right_sf = np.searchsorted(
                dataD["ch5"][1], dataD["ch0"][1][0] + atomDic[atom][1]
            )
            # index of the first syncFastScan that we consider
            left_sf = np.searchsorted(
                dataD["ch5"][1], dataD["ch0"][1][0] + atomDic[atom][0]
            )

            # time after which we start to consider data
            timeInit = dataD["ch5"][1][left_sf] - 1e-9
            # time after which we stop to consider data
            timeEnd = dataD["ch5"][1][right_sf] - 1e-9

            for ch in dataVD.keys():
                # index of the first time after or equal to timeInit
                left = np.searchsorted(dataD[ch][1], timeInit)
                # index of the time before timeEnd
                right = np.searchsorted(dataD[ch][1], timeEnd)

                dataVD[ch][1] = np.append(
                    dataVD[ch][1], dataD[ch][1][left:right])

        return dataVD

    def atomDicCleaner(self, atomDic: dict, specDuration):
        for key in atomDic:
            specNr = (atomDic[key][1] - atomDic[key][0]) // (specDuration)
            atomDic[key][1] = specNr * specDuration + atomDic[key][0]
        return atomDic


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
        kcCounts: int = 4000,
        no=1,
    ):
        print("Post Selecting " + filename + filetype)

        coolingTime = 25e-3
        fsdelay = 0.7e-6
        photonGate = [0, coolingTime]

        # ------ We get the data ------
        file = path + filename + filetype
        filedata = h5py.File(file, mode="r")
        atomnumber = int(len(filedata) / 8)

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
        print(atomList)

        for i in tqdm(atomList, file=sys.stdout):
            # Data loading
            (dataDic) = self.load.data_loading(path, filename, i)

            # We calculate the number of photons per trial during the cooling, optical pumping and state detection periods for the Short cavity

            """
                Really we just count all the all the counts in the short cavity for a run, and then we save:
                    counts in KC -> dataPhotonKC
                    times for an atom in KC -> dataTimeKC

                timeStamps are all the time stamps in seconds
            """
            # Can prolly just leave out the slicing, but let's just keep it here just in case
            # timeStamps = dataDic[self.syncFast][1]

            # --- Short Cavity --- #

            dataPhotonKC = []
            dataTimeKC = []
            timeStamps = dataDic[self.syncFast][1][1:-1]
            for timeStamp in timeStamps:
                """
                    Idea here is that all the timeStamps for events, be it detector events from the
                    short or long cavity, pulses we send to the qTau via the FPGA are being stored as
                    time stamps, with of course a certain delay in between the events.

                    Counting up how many timeStamps we've had in an inteval, will tell us how
                    many events there were.

                """
                timeStamp = timeStamp + fsdelay

                left1, right1 = np.searchsorted(
                    dataDic[self.kcH][1], [timeStamp, timeStamp + coolingTime]
                )
                left2, right2 = np.searchsorted(
                    dataDic[self.kcV][1], [timeStamp, timeStamp + coolingTime]
                )

                # save total number of counts in both cavities at a certain time stamp
                dataPhotonKC.append(right1 - left1 + right2 - left2)
                # saves time stamp
                dataTimeKC.append(timeStamp)

            current_dataPhoton_grouped = []
            current_dataTime_grouped = []
            for current in range(0, len(dataPhotonKC) - no, no):
                """
                    Groups up KC counts
                """
                current_dataPhoton_grouped.append(
                    sum(dataPhotonKC[current: current + no]) / no
                )

            dataPhoton_grouped = dataPhoton_grouped + current_dataPhoton_grouped

            for current in range(0, len(dataPhotonKC) - no, no):
                """
                    Groups up KC timestamps
                """
                current_dataTime_grouped.append(dataTimeKC[current])
            dataTime_grouped = dataTime_grouped + current_dataTime_grouped

            # --- Long Cavity --- #

            dataPhotonLC = []
            dataTimeLC = []
            # We calculate the number of photons per trial for the Long Cavity
            for k, timeStamp in enumerate(dataDic[self.syncFast][1][:-1]):
                timeStamp = timeStamp + fsdelay
                left1, right1 = np.searchsorted(
                    dataDic[self.lcH][1],
                    [timeStamp + photonGate[0], timeStamp + photonGate[1]],
                )
                left2, right2 = np.searchsorted(
                    dataDic[self.lcV][1],
                    [timeStamp + photonGate[0], timeStamp + photonGate[1]],
                )
                dataPhotonLC.append(right1 - left1 + right2 - left2)
                dataTimeLC.append(timeStamp)

            current_dataPhoton_groupedLC = [
                sum(dataPhotonLC[current: current + no]) / no
                for current in range(0, len(dataPhotonLC) - no, no)
            ]
            dataPhoton_groupedLC = dataPhoton_groupedLC + current_dataPhoton_groupedLC
            current_dataTime_groupedLC = [
                dataTimeLC[current] for current in range(0, len(dataPhotonLC) - no, no)
            ]
            dataTime_groupedLC = dataTime_groupedLC + current_dataTime_groupedLC

            # print("Mean ph number = ", np.mean(current_dataPhoton_grouped))

            wt_kc = 0.6 * kcCounts  # wt = witness threshold
            wt_lc = -1  # wt = witness threshold
            # twot = 1.4 * kcCounts  # twot = two atom threshold
            twot = 2 * kcCounts  # twot = two atom threshold
            inAtom = False
            atomIn_index = 0
            atomOut_index = 0

            """
                Now we want to figure out if the counts are in the right range
                so what we do is we're going to iterate over all the grouped
                counts and
            """
            for n, j in enumerate(current_dataPhoton_grouped[1:]):
                n = n + 1
                # check if our counts are above the witness threshold and below the two atom threshold
                if inAtom == False and j >= wt_kc and j >= wt_lc and j <= twot:
                    """
                        if that is the case then we have an atom and detect that
                        an atom has entered the cavity, and we thus set atomIn_index to the index
                        where our atom is in and set inAtom to True since
                        there is an atom in the cavity
                    """
                    atomIn_index = n
                    inAtom = True

                if inAtom:
                    """
                        if our atom is still in, then we just update our atomOut
                        index to the current iteration
                    """
                    atomOut_index = n

                if inAtom == True and (j < wt_kc or j < wt_lc):
                    """
                        if our atom for some reason though gets below the witness threshold
                        we exit out of our loop and and tell the code that we lost the atom
                        by setting the atomOut index to our latest iteration
                    """
                    atomOut_index = n
                    inAtom == False
                    break

                if inAtom == True and (j > twot):
                    """
                        Something similar happens for the case where are atom (or maybe there were two)
                        goes/go above the two atom threshold, which leads to an early exit aswell
                    """
                    atomOut_index = atomIn_index
                    inAtom == False
                    break

            atomInNog2_index = (
                atomIn_index  # as if there would have been no g2 atom selector
            )
            atomOutNog2_index = atomOut_index

            try:
                """
                    now we need to put allat in something we can iterate over in
                    order to be able to plot it later

                    logic here is:
                        timeStamps[atomIn_index] -> gives a time
                        that time - start time from syncSlow signal
                        => gives time when atom entered the cavity


                        timeStamps[atomOut_index] -> gives a time
                        that time - start time from syncSlow signal
                        => gives total lifetime of atom
                """
                # --- Normal Times --- #
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

                # --- Shit dawg ion know --- #
                atomInNog2.append(
                    current_dataTime_grouped[atomInNog2_index]
                    - dataDic[self.syncSlow][1][0]
                )
                atomInNog2Histo.append(
                    current_dataTime_grouped[atomInNog2_index])

                atomOutNog2.append(
                    current_dataTime_grouped[atomOutNog2_index]
                    - dataDic[self.syncSlow][1][0]
                )
                atomOutNog2Histo.append(
                    current_dataTime_grouped[atomOutNog2_index])
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
        """
            Only select the ones where the duration inside of the cavity is
            above a certain threshold
        """
        goodAtomsDF = atomDF[(atomDF["atomsDuration"] >= self.adt)]
        goodAtomsDic = {
            i: [goodAtomsDF["atomsIn"][i], goodAtomsDF["atomsOut"][i]]
            for i in list(goodAtomsDF.index)
        }

        # The conditions for good atoms selection are saved in a data frame
        condsDF = pd.DataFrame()
        condsDF["Conditions"] = ["Single atom time threshold (s)"]
        condsDF["Bounds"] = [self.adt]

        # %% ------ We plot the data ------
        plt.close("all")

        f = plt.figure(filename, figsize=[17, 14])
        f.suptitle("%s, atom %d, binning = %d" % (filename, i, no))

        ax1 = f.add_subplot(211)
        ax2 = f.add_subplot(212)

        # --- kc counts plot --- #
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
        ax1.hlines([twot], atomInHisto[0], atomOutHisto[-1],
                   color="tab:red", alpha=0.2)

        # --- lc counts plot --- #
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
            """
                if an atom is in the cavity and lives long enough the
                background will be dyed in a color in specific color
            """
            if atomOutHisto[i] - atomInHisto[i] >= self.adt:
                ax1.axvspan(
                    atomInHisto[i], atomOutHisto[i], alpha=0.5, color="tab:purple"
                )
                ax2.axvspan(
                    atomInHisto[i], atomOutHisto[i], alpha=0.5, color="tab:purple"
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
            f.savefig(path + filename + ".png")

        return goodAtomsDic, atomInHisto, atomOutHisto

    def getTrapTimes(self, goodAtomsDic, atomInHisto, atomOutHisto):
        list_trappingDuration = []
        for key in goodAtomsDic:
            list_trappingDuration.append(
                goodAtomsDic[key][1] - goodAtomsDic[key][0])

        averageTrapTime = np.mean(list_trappingDuration)
        averageTrapTime_err = np.std(list_trappingDuration) / np.sqrt(
            np.size(list_trappingDuration)
        )

        print(
            "Average single atom trapping time: (%.2f +/- %.2f)s"
            % (averageTrapTime, averageTrapTime_err)
        )

        trappingProbability = len(
            list_trappingDuration) / len(atomInHisto) * 100

        print("Atom trapping probability : %d %%" % (trappingProbability))

        dutyCycle = (
            sum(list_trappingDuration) /
            (atomOutHisto[-1] - atomInHisto[0]) * 100
        )

        print("Duty cycle: %d %%" % (dutyCycle))
        return averageTrapTime, averageTrapTime_err, trappingProbability, dutyCycle


if __name__ == "__main__":
    """
        This script analyzes the trap time and trapping probability
    """

    # --- Definition of plotting stuff --- #
    mot_load_times = [0.6, 0.7, 0.8, 0.9, 1.0]
    trap_times = []
    trap_times_stderr = []
    trap_probabilities = []

    # --- Definition of data sources and destinations --- #
    path = "./data/"
    measurements = "old_trapping"

    if type(measurements) is not list:
        measurements = [measurements]

    filetype: str = ".h5"

    # --- Begin Analysis --- #

    analysis = AtomAnalysis()

    # --- Data evaluation --- #
    for measurement in measurements:
        goodAtomsDic, atomInHisto, atomOutHisto = analysis.dataEv_postSelection(
            path, measurement, filetype
        )
        trap_time, trap_time_stderr, trap_probability, _ = analysis.getTrapTimes(
            goodAtomsDic, atomInHisto, atomOutHisto
        )
        trap_times.append(trap_time)
        trap_times_stderr.append(trap_time_stderr)
        trap_probabilities.append(trap_probability)

    trap_times = np.nan_to_num(trap_times, nan=0.0)
    trap_times_stderr = np.nan_to_num(trap_times_stderr, nan=0.0)

    plt.close("all")
    plt.style.use("seaborn-v0_8-whitegrid")
    # Create the figure and subplots
    fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True, figsize=(8, 6))
    fig.suptitle("MOT Load Time Effects", fontsize=14, fontweight="bold")
    # Top subplot: Trap Times with error bars
    ax0.errorbar(
        mot_load_times,
        trap_times,
        yerr=trap_times_stderr,
        fmt="o-",
        capsize=4,
        color="tab:blue",
        ecolor="gray",
    )
    ax0.set_title("Trap Times", fontsize=12)
    ax0.set_ylabel("Trap Time (s)")
    ax0.grid(True)
    # Bottom subplot: Trap Probability
    ax1.plot(
        mot_load_times, trap_probabilities, marker="o", linestyle="-", color="tab:green"
    )
    ax1.set_title("Trap Probability", fontsize=12)
    ax1.set_xlabel("MOT Load Time (s)")
    ax1.set_ylabel("Probability (a.u.)")
    ax1.grid(True)
    # Tight layout with space for suptitle
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig("new_trapping.png")
