---
id: "Laser_Beam_Preparation"
date: "2025-05-16"
subtitle: 
tags:
    - light
    - experiment
---

# Laser_Beam_Preparation

Light which is applied to the atom mnust be adjusted in at least four degrees of freedom:

    1. The Frequency
    2. The Polarization
    3. The optical Power
    4. The turn on and off time

These properties are set in our experiment via two beam paths that are repeatedly built.

The first one superimposes the laser light with a reference light in order to stabilize the laser light frequency.
The second path uses [AOMs](general_knowledge/AOM.md) in a double-pass configuration to fine tune the laser frequency, optical power, and to switch the light on and off on a sub-microsecond timescale,
after which it gets coupled into a [PM fiber](general_knowledge/Fibers.md) which guides the light to the optical table with the cavities.

The seperation of the two optical tables is intended, as the heat dissipation by all the components, may cause thermal drifts. Furthermore it makes maintanance work on the laser lights easier, as working
on the experimental table is quite difficult due to the spatial constrains, and they increase the risk of causing damage to the experimental setup.

## Frequency Lock of the Laser Sources

Problem with `ECDLs` or `External Cavity Diode Lasers` is that they tend to spectrally drift over time over a spectral range far broader than than the natural linewidth of the employed $^{87}Rb$ transitions.
This then means that we have to lock the lasers to some sort of frequency reference, which in our case is done via an optical `frequency comb`.

As the name suggests, the spectral structure of a frequency comb is that of, well a comb. The combs are spaced $250 MHz$ apart with an offest frequency of $-20 MHz$ over a range of $(760 - 805) nm$.
Therefore the aboslute positions of the comb teeth are at $f = (-20 + N\cdot250)MHz$ where N is an integer.
![frequency_comb_spectrum.png](assets/imgs/frequency_comb_spectrum.png)

What we do with the frequency comb is, we superimpose it with the laser light, which leads to a beat signal. This signal then gets measured by a photo detector, and then forwarded to the locking electronic with
the corresponding voltage signal. The locking electronic then compares the beat frequency with an adjustable reference radiofrequency and outputs a corresponding control signal.
The signal is then fed back to the ECDL current or grating piezo in order to stabilize the laser frequency.

What would be nice if we could distribute the comb light equally over all lasers that need to be locked in frequency, however with each laser the power per laser becomes small.
What we do instead is we group lasers with a similar target frequency together and filter out the respective comb light spectrum.
![laser_lock_via_frequency_comb_setup.png](assets/imgs/laser_lock_via_frequency_comb_setup.png)

## AOM Track

As the name suggests, `AOMs` are used in this track in a double-pass configuration.
Moreover, the efficiency of the AOM diffraction is a function of the applied RF signal power, which provides an excellent mean for the fine adjustment of the optical power.
Since the diffraction of the laser light is due to the vibrating crystal lattice, the interaction is basically between photons and phonons, meaning the diffracted laser light either lost or gained energy. FOr the first order the change is equivalent ot the applied RF times the Planck constant.

But the AOM can also be used to adjust the laser frequency. 

The setup can be seen here:
![aom_track_setup.png](assets/imgs/aom_track_setup.png)

The way its done here shall be explained now. Light is partially coupled into an AOM track via a $\lambda/2$ waveplate and a `PBS` (`Polarizing Beam Splitter`). This already allows for some `rough` power alignment.
Light is then focused onto an AOM. 

---

`NOTE: The focused beam reduces clipping losses at the AOM aperture and decreases time taken to turn laser pulse on and off. However a smaller beam diameter causes a lower diffraction efficiency. Ergo there's a tradeoff
      between switching time and diffraction efficiency.`

---

After the AOM, the light passes through a $\lambda/4$ waveplate and then an iris which absorbs the 0th diffraction order, so only the 1st one (or more specifically here -1st) is transmitted.
The light then gets reflected, so iris, then $\lambda/4$ waveplate, and AOM, where the direction of diffraction order -1 equals the reverse direction of the initial incoming beam.
The double-pass AOM conﬁguration enables a double energy shift of the laser beam. However, the AOM diffraction ineffciency enters the overall AOM track efficiency quadratically.

The light then gets coupled out of the track by a PBS, as it now has a perpendicular polarization compared to the one it entered the path with, due to the double passage of the $\lambda/4$ waveplate.
The linearly polarized light is then coupled into a PM fiber which guides the light to the cavity table.

---

`NOTE: It is not recommended to couple the laser beam out of the AOM track via transmission through the initial PBS. This is because an imperfect field polarization may lead to a residual beam reflection as opposed to full transmission. Hence, the outcoupled light would be guided towards the previous AOM track, leading to cross talk between two adjacent AOM tracks.`

---

## Laser beam frequencies and the $D_2$-line of $^{87}Rb$

Two ECDLs are used to couple both hyperfine ground states, $F=1$ and $F=2$ to the $5^2P_{3/2}$ manifold. You could use only one ECDL, split the beam and spectrally shift one of the beams with an [EOM](general_knowledge/EOM.md), but this would result in significant optical losses due to the finite efficiency of EOMs. 

In this setup, two lasers are blue-detuned from the unshifted D2-line of rubidium-87 by 295 MHz and 218 MHz. One laser (ECDL TA pro) is locked to a frequency comb at 384.22848 THz with a red detuning of 70 MHz, while the other (ECDL DL pro) is referenced to 384.23523 THz with a red detuning of 62 MHz. 
These lasers are frequency-shifted using double-pass AOM (acousto-optic modulator) tracks, which red-shift the laser light by twice the RF frequency. 
The AOMs have center frequencies and bandwidths of a few hundred MHz, and adjusting the RF signal within this range typically requires re-aligning the AOM for optimal diffraction efficiency.
![laser_beam_frequencies_setup.png](assets/imgs/laser_beam_frequencies_setup.png)


##### References
Niemietz, D. Nondestructive detection of photonic qubits with single atoms in crossed ﬁber cavities. (Technische Universität München,2021), https://mediatum.ub.tum.de/1601857
