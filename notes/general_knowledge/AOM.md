---
id: "AOM"
date: "2025-05-09"
subtitle: 
tags:
    - device
---

# AOM

An `AOM` or `acousto-optic modulator`, also called a Bragg cell, is a device that modulates certain properties of light that pass through a quartz crystal by using sound waves.
An AOM is not to be confused with an [EOM](general_knowledge/EOM.md) as that device relies on the use of electrical fields that get applied to very fancy crystals in order to modulate properties of light fields.

---

The basic setup of an AOM is as follows:
![AOM.png](assets/imgs/AOM.png)

It consists of a piece of glass (usually Quartz) that is stuck between a piezo-electric transducer and an acoustic absober. The `Piezo` is attached so that when we send an electrical signal to it, the transducer vibrates, which creates sound waves in the material. The acoustic waves then create expansion and compression in the glass crystal that change the index of refraction. The grating constant `d` which is then usually observed in `Bragg Reflection` is just the distance between two peaks in the refractive index pattern, which is one full wavelength $\Lambda$ of the acoustic wave: $\Rightarrow d=\Lambda$. Thus we get the following Bragg condition:
$2\Lambda \sin(\theta) = m\lambda$

Typically one uses the first order diffraction as it is the strongest one.
Finally one can use this device to modulate the frequency, phase, intensity and polarization of our light wave.

Most notably the frequency is modulated since the light scatters off of moving planes, and thus leading to a doppler shift of our wave by an equal amount to the acoustic frequency:
$f\rightarrow f + m\Lambda$

This frequency shift can also be understood by the fact that energy and momentum are conserved in the scattering process.
