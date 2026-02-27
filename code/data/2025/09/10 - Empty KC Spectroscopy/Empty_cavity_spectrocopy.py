#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 19 12:11:16 2021

@author: tfrank
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import lmfit

# %% Input variables
scanRange = 5.5  # (GHz) measured with Wavemeter at low scanning frequency
FSR = 1902.4  # Measured with two lasers at different frequ. being in resonance with cavity

conversion_flag = False  # set to False to determine time-frequency-conversion factor, set to True once conversion factor is determined and placed in the conversion section of this script file
saveFit_Flag = True
saveCalibration_Flag = False


# %% FUNCTIONS
# ------ Functions related to curve fitting ------
# Cavity transmission expressions taken from M. Brekenfeld  Thesis (eq. 2.5)
def transmission_cavity(freq, amp1, freqC1, k1, amp2, freqC2, k2):
    t1 = amp1 / (1j * (freq - freqC1) + k1)
    t2 = amp2 / (1j * (freq - freqC2) + k2)
    return abs(t1) ** 2 + abs(t2) ** 2


def transmission_cavity_mod(freq, amp1, freqC1, k1, amp2, freqC2, k2, amp3, freqC3, k3):
    t1 = amp1 / (1j * (freq - freqC1) + k1)
    t2 = amp2 / (1j * (freq - freqC2) + k2)
    t3 = amp3 / (1j * (freq - freqC3) + k3)
    return abs(t1) ** 2 + abs(t2) ** 2 + abs(t3) ** 2


def transmission_cavity_mod5(
    freq,
    amp1,
    freqC1,
    k1,
    amp2,
    freqC2,
    k2,
    amp3,
    freqC3,
    k3,
    amp4,
    freqC4,
    k4,
    amp5,
    freqC5,
    k5,
):
    t1 = amp1 / (1j * (freq - freqC1) + k1)
    t2 = amp2 / (1j * (freq - freqC2) + k2)
    t3 = amp3 / (1j * (freq - freqC3) + k3)
    t4 = amp4 / (1j * (freq - freqC4) + k4)
    t5 = amp5 / (1j * (freq - freqC5) + k5)
    return abs(t1) ** 2 + abs(t2) ** 2 + abs(t3) ** 2 + abs(t4) ** 2 + abs(t5) ** 2


def transmission_cavity_single(freq, amp1, freqC1, k1):
    t1 = amp1 / (1j * (freq - freqC1) + k1)
    return abs(t1) ** 2


def reflection_cavity(freq, freqC, MM_rf, MM_fc, k, k_oc, phi_rf, phi_fc):
    r = MM_rf * np.exp(1j * phi_rf) - (MM_fc * np.exp(1j * phi_fc)) ** 2 * 2 * k_oc / (
        1j * (freq - freqC) + k
    )
    return abs(r / MM_rf * np.exp(1j * phi_rf)) ** 2


# ------ Function to fit the locked cavity reflection and transmission spectrum ------
def dataFit_cavityTR(freqs, tCav, fitMax, fitMin, conversion_flag, ModeSplit):
    # CAVITY TRANSMISSION
    if conversion_flag == True:
        modCavT = lmfit.Model(transmission_cavity)  # model used for fitting
    else:
        modCavT = lmfit.Model(
            transmission_cavity_mod
        )  # model used for fitting conversion/calibration data
    params = modCavT.make_params()

    # Fitting parameter limits definition
    def pset(name, value, min=None, max=None, vary=True):
        params[name].set(value=value, min=min, max=max, vary=vary)

    fitCenter = (fitMax + fitMin) / 2  # estimate location of peaks
    fitrange = fitMax - fitMin

    if conversion_flag == True:
        # Parameters for the cavity locked fit pset(label,start,lowerLimit,upperLimit, fittingBoolean)
        pset("amp1", 0.12, 0, 0.25, True)  # arbitrary amplitude
        pset(
            "freqC1",
            fitCenter - ModeSplit / 2,
            fitCenter - ModeSplit / 2 - 0.05,
            fitCenter - ModeSplit / 2 + 0.05,
            True,
        )  # cavity resonance frequency
        pset("k1", 0.066, 0, 1, True)  # total cavity decay rate
        pset("amp2", 0.02, 0.005, 0.03, True)  # arbitrary amplitude
        pset(
            "freqC2",
            fitCenter + ModeSplit / 2,
            fitCenter + ModeSplit / 2 - 0.05,
            fitCenter + ModeSplit / 2 + 0.05,
            True,
        )  # cavity resonance frequency
        pset("k2", 0.068, 0.060, 0.075, True)  # total cavity decay rate
    else:
        # Parameters for the cavity locked fit pset(label,start,lowerLimit,upperLimit, fittingBoolean)
        pset("amp1", 0.000014, 0, 0.01, True)  # arbitrary amplitude
        pset(
            "freqC1", fitCenter - fitrange / 6, fitMin, fitMax, True
        )  # cavity resonance frequency
        pset("k1", 0.000035, 0, 0.002, True)  # total cavity decay rate
        pset("amp2", 0.000014, 0, 0.01, True)  # arbitrary amplitude
        pset(
            "freqC2", fitCenter + fitrange / 6, fitMin, fitMax, True
        )  # cavity resonance frequency
        pset("k2", 0.000035, 0, 0.002, True)  # total cavity decay rate
        pset("amp3", 0.00006, 0, 0.25, True)  # arbitrary amplitude
        pset("freqC3", fitCenter, fitMin, fitMax, True)  # cavity resonance frequency
        pset("k3", 0.000035, 0, 1, True)  # total cavity decay rate

    res_emptyCavityT = modCavT.fit(tCav, params, freq=freqs, options={"maxfev": 1})
    print(res_emptyCavityT.fit_report())
    return res_emptyCavityT


def dataFit_cavityTR5(freqs, tCav, fitMax, fitMin, conversion_flag, ModeSplit):
    modCavT = lmfit.Model(
        transmission_cavity_mod5
    )  # model used for fitting conversion/calibration data
    params = modCavT.make_params()

    # Fitting parameter limits definition
    def pset(name, value, min=None, max=None, vary=True):
        params[name].set(value=value, min=min, max=max, vary=vary)

    fitCenter = (fitMax + fitMin) / 2  # estimate location of peaks
    fitrange = fitMax - fitMin

    # Parameters for the cavity locked fit pset(label,start,lowerLimit,upperLimit, fittingBoolean)
    pset("amp1", 0.0007, 0, 0.01, True)  # arbitrary amplitude
    pset("freqC1", fitCenter, fitMin, fitMax, True)  # cavity resonance frequency
    pset("k1", ModeSplit / 4.17, 0, 0.0002, True)  # total cavity decay rate

    pset("amp2", 0.0007, 0, 0.01, True)  # arbitrary amplitude
    pset(
        "freqC2", fitCenter + ModeSplit / 2, fitMin, fitMax, True
    )  # cavity resonance frequency
    pset("k2", ModeSplit / 4.17, 0, 0.0002, True)  # total cavity decay rate

    pset("amp3", 0.00016, 0, 0.25, True)  # arbitrary amplitude
    pset(
        "freqC3", fitCenter - ModeSplit / 2, fitMin, fitMax, True
    )  # cavity resonance frequency
    pset("k3", ModeSplit / 4.17, 0, 0.0002, True)  # total cavity decay rate

    pset("amp4", 0.00016, 0, 0.25, True)  # arbitrary amplitude
    pset(
        "freqC4", fitCenter + ModeSplit, fitMin, fitMax, True
    )  # cavity resonance frequency
    pset("k4", ModeSplit / 4.17, 0, 0.0002, True)  # total cavity decay rate

    pset("amp5", 0.00016, 0, 0.25, True)  # arbitrary amplitude
    pset(
        "freqC5", fitCenter - ModeSplit, fitMin, fitMax, True
    )  # cavity resonance frequency
    pset("k5", ModeSplit / 4.17, 0, 0.0002, True)  # total cavity decay rate

    res_emptyCavityT = modCavT.fit(tCav, params, freq=freqs, options={"maxfev": 1})
    print(res_emptyCavityT.fit_report())
    return res_emptyCavityT


def dataFit_cavityRF(freqs, rCav, fitMax, fitMin, ModeSplit):
    # CAVITY Rreflection
    modCavR = lmfit.Model(reflection_cavity)  # model used for fitting
    params = modCavR.make_params()

    # Fitting parameter limits definition
    def pset(name, value, min=None, max=None, vary=True):
        params[name].set(value=value, min=min, max=max, vary=vary)

    fitCenter = (fitMax + fitMin) / 2  # estimate location of peaks
    fitrange = fitMax - fitMin

    # Parameters for the cavity locked fit pset(label,start,lowerLimit,upperLimit, fittingBoolean)
    pset("MM_rf", 0.99, 0, 1, True)  # arbitrary amplitude
    pset("MM_fc", 0.882, 0.85, 0.95, True)  # cavity resonance frequency
    pset("k_oc", 116 / 2 * 0.85, 0, 70, False)  # total cavity decay rate
    pset("freqC", 0, -10, +10, True)  # cavity resonance frequency
    # pset("freqC", fitCenter-ModeSplit/2, fitCenter-ModeSplit/2-0.05, fitCenter-ModeSplit/2+0.05, True)  # cavity resonance frequency
    pset("k", 115 / 2, 0, 70, False)  # total cavity decay rate
    pset("phi_rf", 0.00, 0, 2 * np.pi, False)  # total cavity decay rate
    pset("phi_fc", -0.05, -np.pi, 2 * np.pi, True)  # total cavity decay rate
    # pset("amp3", 0.14,0, 0.25, True)  # arbitrary amplitude
    # pset("freqC3", fitCenter, fitMin, fitMax, True)  # cavity resonance frequency
    # pset("k3", 0.0016, 0, 1, True)  # total cavity decay rate

    res_emptyCavityR = modCavR.fit(rCav, params, freq=freqs, options={"maxfev": 1})
    print(res_emptyCavityR.fit_report())
    return res_emptyCavityR


def get_max_cavityT(Npeak, cavityT):
    I = len(cavityT["Ampl"]) // Npeak
    cavityTpeaks = np.array(cavityT["Ampl"][0 : I * Npeak]).reshape((Npeak, I))
    list_cavityTmax = []

    for i in range(Npeak):
        maxIndex = np.argmax(cavityTpeaks[i])
        list_cavityTmax.append(cavityT["Freq"][maxIndex + I * i])
    return list_cavityTmax


def get_min_cavityR(Npeak, cavityR):
    I = len(cavityR["Ampl"]) // Npeak
    cavityRpeaks = np.array(cavityR["Ampl"][0 : I * Npeak]).reshape((Npeak, I))
    list_cavityRmin = []

    for i in range(Npeak):
        minIndex = np.argmin(cavityRpeaks[i])
        list_cavityRmin.append(cavityR["Freq"][minIndex + I * i])
    return list_cavityRmin


def saveFitData(directory, fit_emptyCavityT, fit_emptyCavityR):
    if saveFit_Flag == True:
        if "linewidth_2.txt" in os.listdir(directory):
            with open(directory + "linewidth_2.txt", "a+") as file:
                file.write(
                    "{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}, {9}".format(
                        fit_emptyCavityT.best_values["freqC1"],
                        fit_emptyCavityT.best_values["amp1"],
                        fit_emptyCavityT.best_values["k1"],
                        Finesse1,
                        fit_emptyCavityT.best_values["freqC2"],
                        fit_emptyCavityT.best_values["amp2"],
                        fit_emptyCavityT.best_values["k2"],
                        Finesse2,
                        modeSplit,
                        fit_emptyCavityR.best_values["MM_fc"],
                    )
                )
                file.write("\n")
            file.close
        else:
            with open(directory + "linewidth_2.txt", "w+") as file:
                file.write("fit data from" + fileName_cavityT)
                file.write("\n")
                file.write(
                    "resonance position 1 (GHz), lw1 (GHz), Finesse1, resonance position 2 (GHz), lw2 (GHz), Finesse 2, splitting"
                )
                file.write("\n")
                file.write(
                    "{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}, {9}".format(
                        fit_emptyCavityT.best_values["freqC1"],
                        fit_emptyCavityT.best_values["amp1"],
                        fit_emptyCavityT.best_values["k1"],
                        Finesse1,
                        fit_emptyCavityT.best_values["freqC2"],
                        fit_emptyCavityT.best_values["amp2"],
                        fit_emptyCavityT.best_values["k2"],
                        Finesse2,
                        modeSplit,
                        fit_emptyCavityR.best_values["MM_fc"],
                    )
                )
                file.write("\n")
            file.close()


def saveCalibrationData(directory, conversion):
    if saveCalibration_Flag == True:
        if "conversion.txt" in os.listdir(directory):
            with open(directory + "conversion.txt", "a+") as file:
                file.write("{0}".format(conversion))
                file.write("\n")
            file.close
        else:
            with open(directory + "conversion.txt", "w+") as file:
                file.write("data for time to frequency conversion")
                file.write("\n")
                file.write("conversion factor")
                file.write("\n")
                file.write("{0}".format(conversion))
                file.write("\n")
            file.close()


# %%LOAD DATA
directory = "Z:/Results/2025/09/10 - Empty KC Spectroscopy/"
fileName_cavityT = "kc_transmission_300MHzModulation.csv"
fileName_cavityR = "C2--kc_reflection_pi_mode--00001.csv"  # contains reflection signal
# fileName_cavityRBG = 'C1--reflectionoffset--00000.csv' #contains1Background signal
fileName_cavityTBG = (
    "kc_transmission_300MHzModulation.csv"  # contains Background signal
)
data_cavityT = pd.read_csv(fileName_cavityT, delimiter=",", skiprows=4)
data_cavityR = pd.read_csv(fileName_cavityR, delimiter=",", skiprows=4)
# data_cavityRBG=pd.read_csv(directory+fileName_cavityRBG, delimiter=',', skiprows=4)
data_cavityTBG = pd.read_csv(fileName_cavityTBG, delimiter=",", skiprows=4)

# %%Convert time to frequency using eom modulated side bands


# rename time to frequency
data_cavityT.columns = ["Freq", "Ampl"]
data_cavityR.columns = ["Freq", "Ampl"]

# Normalization
data_cavityR["Ampl"] = data_cavityR["Ampl"]
data_cavityR["Ampl"] = data_cavityR["Ampl"] / (np.mean(data_cavityR["Ampl"][0:1000]))
data_cavityR["Freq"] = data_cavityR["Freq"]  # + 0.0022
data_cavityT["Ampl"] = data_cavityT["Ampl"] - np.mean(data_cavityT["Ampl"][0:1000])


plt.plot(data_cavityR["Freq"], data_cavityR["Ampl"], rasterized=True)
plt.plot(data_cavityT["Freq"], data_cavityT["Ampl"], rasterized=True)
plt.show()
# %% DATA FITTING
list_linewidth = []
list_splitting = []
list_MM_fc = []
list_MM_fc_phi = []
list_MM_fr = []
splitting_flag = True
if splitting_flag == True:
    list_conversion1 = []
    list_conversion2 = []
    list_cavityTmax = get_max_cavityT(50, data_cavityT)
    list_cavityTmax2 = []
    list_cavityRmin = get_min_cavityR(50, data_cavityR)
    plt.rcParams.update({"font.size": 14})
    spec_fig = plt.figure(figsize=[12, 8])
    spec_fig.suptitle("Empty Cavity Spectroscopy")
    ax1 = spec_fig.add_subplot(2, 1, 1)
    ax2 = spec_fig.add_subplot(2, 1, 2)
    print(list_cavityTmax)
    # Set fit interval symetrically around the resonance to be fit
    i = 0
    for freq in list_cavityTmax[0:50]:
        # print(freq)
        # Set fit interval symetrically around the resonance to be fit

        ModeSplit = 0.001  # 00005/0008/00011
        fitC = freq
        fitI = ModeSplit

        fitMin = fitC - fitI / 2
        fitMax = fitC + fitI / 2
        dfFit = data_cavityT[(data_cavityT["Freq"] > fitMin)]
        dfFit = dfFit[(dfFit["Freq"] < fitMax)]
        fit_emptyCavityT = dataFit_cavityTR(
            dfFit["Freq"], dfFit["Ampl"], fitMax, fitMin, conversion_flag, ModeSplit
        )

        conversion = abs(
            fit_emptyCavityT.best_values["freqC1"]
            - fit_emptyCavityT.best_values["freqC2"]
        )  # *1000
        list_conversion1.append(conversion)
        # conversion = ModeSplit/2
        # print(conversion)

        linewidth = fit_emptyCavityT.best_values["k3"] * 2 / conversion * 600
        list_linewidth.append(linewidth)

        dfFitR = data_cavityR[(data_cavityR["Freq"] > fitMin)]
        dfFitR = dfFitR[(dfFitR["Freq"] < fitMax)]

        if np.mod(i, 2) == 0:
            fitCR = list_cavityRmin[i]
            fit_emptyCavityR = dataFit_cavityRF(
                (dfFitR["Freq"] - fitCR) / conversion * 600,
                dfFitR["Ampl"],
                fitMax,
                fitMin,
                ModeSplit,
            )
            if np.mod(i, 5) == 0:
                ax2.plot(
                    (dfFitR["Freq"] - fitCR) / conversion * 600,
                    dfFitR["Ampl"],
                    rasterized=True,
                )

            modCavR = lmfit.Model(reflection_cavity)
            if np.mod(i, 5) == 0:
                ax2.plot(
                    (dfFitR["Freq"] - fitCR) / conversion * 600,
                    modCavR.eval(
                        freq=(dfFitR["Freq"] - fitCR) / conversion * 600,
                        params=fit_emptyCavityR.params,
                        rasterized=True,
                    ),
                    color="r",
                    zorder=2,
                )
            list_MM_fc.append(fit_emptyCavityR.best_values["MM_fc"])
            list_MM_fc_phi.append(fit_emptyCavityR.best_values["phi_fc"])
            list_MM_fr.append(fit_emptyCavityR.best_values["MM_rf"])

        if np.mod(i, 5) == 0:
            ax1.plot((dfFit["Freq"] - fitC) / conversion * 600, dfFit["Ampl"])
        modCavT = lmfit.Model(transmission_cavity_mod)

        ax1.plot(
            (dfFit["Freq"] - fitC) / conversion * 600,
            modCavT.eval(freq=dfFit["Freq"], params=fit_emptyCavityT.params),
            color="r",
            zorder=2,
        )

        # plt.plot((dfFitR['Freq']-fitC)/conversion*500,
        #         reflection_cavity((dfFitR['Freq']-fitC)/conversion*500,0,0.99,0.93,68/2,58/2,0.05,0))

        ax1.set_ylabel("Intensity (a.u.)")
        ax2.set_xlabel("Frequency (MHz)")
        ax2.set_ylabel("Intensity (a.u.)")

        i = i + 1
    splitting_flag = False

    print("###########splittingFitend###########")


# %%
# if conversion_flag == False:
#     print(list_conversion)
#     print(np.average(list_conversion))


# print(list_splitting)
# mean_splitting = np.average(list_splitting)
# std_of_mean_splitting = np.std(list_splitting)/np.sqrt(len(list_splitting))
# print(mean_splitting)
# print(std_of_mean_splitting)
# print(np.std(list_splitting))


# print(list_linewidth)
mean_linewidth = np.average(list_linewidth)
std_of_mean_linewidth = np.std(list_linewidth) / np.sqrt(len(list_linewidth))
print("Linewidth = %.2f(%.2f) MHz" % (mean_linewidth, std_of_mean_linewidth))
# print(std_of_mean_linewidth)
# print(np.std(list_linewidth))
ax1.set_title(
    # "KC transmission with EOM sidebands @ 300 MHz\n"
    "Transmission through Cavity \n"
    + r"Linewidth FWHM; $2\kappa$= %.2f(%.f) MHz"
    % (mean_linewidth, std_of_mean_linewidth * 100)
)

# print(list_linewidth)
mean_MM_fc = np.average(list_MM_fc)
mean_MM_fc_phi = np.average(list_MM_fc_phi)
mean_MM_fr = np.average(list_MM_fr)
std_of_mean_MM_fc = np.std(list_MM_fc) / np.sqrt(len(list_MM_fc))
std_of_mean_MM_fc_phi = np.std(list_MM_fc_phi) / np.sqrt(len(list_MM_fc_phi))
std_of_mean_MM_fr = np.std(list_MM_fr) / np.sqrt(len(list_MM_fr))
print(mean_MM_fc, std_of_mean_MM_fc)
print(
    "Mode matiching fiber-cavity = %.3f(%.f) MHz"
    % (mean_MM_fc, std_of_mean_MM_fc * 1000)
)

print(
    "Angle Mode matiching fiber-cavity = %.3f(%.f) MHz"
    % (mean_MM_fc_phi, std_of_mean_MM_fc_phi * 1000)
)

print(
    "Mode matiching fiber-reflected light = %.3f(%.f) MHz"
    % (mean_MM_fr, std_of_mean_MM_fr * 1000)
)
# print(std_of_mean_linewidth)
# print(np.std(list_linewidth))
ax2.set_title(
    "Reflection \n"
    + r"Mode matching fiber-cavity; $\mu_{fc}$ = %.3f(%.f)"
    % (mean_MM_fc, std_of_mean_MM_fc * 1000)
    + "\n"
    + r"Mode matching fiber-reflected light; $\mu_{fr}$ = %.3f(%.f)"
    % (mean_MM_fr, std_of_mean_MM_fr * 1000)
)
ax2.set_ylim(0, 1.1)
plt.tight_layout()
print(f"Reflection for |0>: {np.min(dfFitR['Ampl'])}")
G0_kc = 2 * np.pi * 0.0386  # 2-1' splitting
Kappa = 2 * np.pi * 0.058
Kappa_oc = Kappa * 0.85
Gamma_5P32_5S = 2 * np.pi * 0.006065 / 2
Coopoerativity = G0_kc**2 / (2 * Kappa * Gamma_5P32_5S)
Mu_rf = 0.88
Mu_fc = 0.88

print(
    f"Phase shift for |0>: {np.arctan2(0, (Kappa * Gamma_5P32_5S**2 * (Kappa - 2 * Kappa_oc)))}"
)
plt.savefig("Normal_mode_spectroscopy.svg")
plt.show()
