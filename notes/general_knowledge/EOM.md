---
id: "EOM"
date: "2025-05-09"
subtitle: 
tags:
    - device
---

# EOM

An `EOM` or `electro-optical modulator` is an optical device that is there in order to modulate certain properties of a beam of light by the application of an electrical field to a crystal, which sits at the heart of this device. 
An EOM is not to be confused with an [AOM](general_knowledge/AOM.md) as that device relies on the use of acoustic waves that get applied to a piece of glass in order to modulate properties of light fields.

Modulation in this case means, that we change a property (e.g. Amplitude, Phase, Frequency) in some sense by a certain amount

---

An EOM looks like this:
![EOM.png](assets/imgs/EOM.png)

It has an input, output where we stick our fibres in for our light to pass through. And then there is the `RF` port, where we can attach one of our cables to in order to pass a signal onto the crystal.

The theory behind the EOM is, that in it sits a crystal, whose refractive index is a function of the local electric field (i.e. our signal we send to the RF port). If we now put our crystal in iside of a parallel plate capacitor, link that to our RF port, we can control the electrical field in the capacitor with our signal, and by that extension the refractive index of our crystal. 

Changing the refractive index of our crystal, then means we change the speed at which light is able to travel through our crystal (as the speed at which light travels through a medium is given by $v = \frac{c}{n}$). But since the phase of the light leaving the crystal is directly proportional to the length of time it takes for light to pass through it, we thus control the phase of the laser light exiting an EOM by changing the electric field in the crystal via our RF port.

##### References
Wikipedia contributors. (2024, April 18). Electro-optic modulator. Wikipedia. https://en.wikipedia.org/wiki/Electro-optic_modulator#EOM_technologies
