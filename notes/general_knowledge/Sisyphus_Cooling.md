---
id: "Sisyphus_Cooling"
date: "2025-05-09"
subtitle: 
tags:
    - experiment
    - cooling
---

# Sisyphus_Cooling

This is going to be a summarization of what `Sisyphus Cooling` is and how it works and whatnot
This summarization is going to be based on the book `Atomic Physics by Christopher J.Foot (2005)`

---

## Method

Assuming that we have an atom that is very `hot`. We can determine how `hot` our atom is by looking at how much it vibrates around a mean position (which in essence is all that temperature is).
Now we're looking to reduce this vibration, which in another sense is just kinetic Energy.
In order to perform Sisyphus Cooling we're sending two counterpropagating light beams with perpendicular polarization in the direction of our atom.
This will create a standing wave in the polarization regime that is alternating between left circular, linear, right circular, and linear polarized light across the wave.
![sisyphus_interference_standing_wave.png](assets/imgs/sisyphus_interference_standing_wave.png)

Now this is important since when we look at our energy level scheme of our atom:
![sisyphus_atom_energy_scheme.png](assets/imgs/sisyphus_atom_energy_scheme.png)

We see that our ground states at $J=\frac{1}{2}$ couple differently strong to our different circular polarized light fields thus giving us different light shifts depending on in which ground state we are:
![sisyphus_light_shift_energy.png](assets/imgs/sisyphus_light_shift_energy.png)

These shifts depend on the `Clebsch-Gordan coefficients`, and need to be calculated.

The idea is now that our atom has a certain amount of kinetic energy, it can take some of that kinetic energy to move up the potential hill, convert that energy into potential energy, and thus slow down a bit. But if we were to do nothing that potential energy would just get re-transferred back into kinetic energy. If we now use a laser that is resonant to the transition when the energy of the atom is on top of the hill, we can excite our atom to an excited state and it will then decay via spontaneous emission to the other ground state and we keep repeating that process, until the atom doesn't have enough kinetic energy to move up that hill anymore. 
![sisyphus_laser_cooling_mechanism.png](assets/imgs/sisyphus_laser_cooling_mechanism.png)

---

### Limits

Sisyphus cooling stops working when the loss in energy in going from teh top of a hill (in the potential energy) to the valley bottom is balanced by the recoil energy acquired in spontaneous emission, $U_0 \simeq E_r$.
For this condition there is no net loss of energy in optical pumping between $M_J$ states. Thus the lowest temperatures achieved are equivalent to a few recoil energy, $T\simeq\frac{E_r}{k_b}$. At this recoild limit the temperature is given by:
![sisyphus_recoil_limit.png](assets/imgs/sisyphus_recoil_limit.png)
