---
id: "Waveplates_for_polarization_basis"
date: "2025-12-19"
subtitle: 
tags:
---

# Waveplates_for_polarization_basis

We got two waveplates to set the basis of our polarization qubit. They're controlled by the ELL14 motorized mounts by thorlabs which are controlled by a daughterbord which connects both waveplates to the main board which then connects to the pc. 

They have the following addresses in the ello software:

- The L/2 waveplate has address 0
- The L/4 waveplate has address 1

In order to control them I have coded a gui, where one needs to first of all find out all the waveplate settings in order to resolve each polarization.
It can be found in the control folder. I also included the settings with which I'm able to resolve all polarizations (state: 19.12.2025). Maybe I'll have to update them, although I don't think that's going to be the case, since I have a polarizer right in front of the WPs and the chance of that thing tweaking are I think very low.

The gui uses this [piece of code](https://github.com/roesel/elliptec) (thank goodness some mf actually coded out a comprehensible piece of software to controls these POSs). In order to start the code make sure you have pyqt5 installed, and just run the main file. 

## Note

Make sure the usb cable plugged in, that you have the right usb-port selected, and also that each waveplate has distinctive addresses. Otherwise the software will not work properly.
