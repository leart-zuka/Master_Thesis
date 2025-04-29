---
id: "quantum_gate_between_flying_photon_and_trapped_atom"
date: "2025-04-29"
subtitle: 
tags:
    - paper
---

# quantum_gate_between_flying_photon_and_trapped_atom

This is going to be my summarization of the paper `A quantum gate between a flying optical photon and a single trapped atom` 
(DOI:[10.1038/nature13177](https://www.nature.com/articles/nature13177))

## Intro

We know that quantum computers and quantum networks are very useful in the way, that quantum computers allows for a large increase in computability
if it can be scaled to large number of qubits, and quantum networks based on photons guarentee unbreakable encryption but it has to be scaled to high rates
and large distances.

One way to get around this issue is to transition to a hybrid system of qubits and matter qubits, in order to solve scalability in each field-that of communication by use of quantum
repeaters, and that of computation by use of an optical interconnect between smaller quantum processors.

In order to link these two systems however it requires the construction of a quantum gate.
This paper presents a quantum gate between the spin state of a single trapped atom and the polarization state of an optical photon contained in a faint laster pulse.

## Idea

The method used in this paper is based on cavity quantum electrodynamics. When a photon interacts with a cavity containing a single resonant emitter, it experiences a phase shift which depends on the coupling strength.
Here we're using a single Rubidium atom, which is trapped at the centre of an overcoupled cavity. With the parameters of that paper, the conditional phase shift induced on the photon is $\pi$, which is the prerequesite for this quantum gate presented in this paper With the parameters of that paper, the conditional phase shift induced on the photon is $\pi$, which is the prerequesite for this quantum gate presented in this paper.

Because looking at the energy scheme of our Rubidium atom inside of a strong electric field, we see the following:
![rubidium_energy_scheme_stark_shift.png](assets/imgs/rubidium_energy_scheme_stark_shift.png)
Thus the coupling is only strong when the atom (a) is in a state $\ket{\uparrow^a}$ and photons (p) of right-circular polarization $\ket{\uparrow^p}$ are reflected (green arrow and sphere in figure). A pi shift in this case
would be the same as a sign change between the atomic and photonic qubit:
![pi_shift_atom_photon_state.png](assets/imgs/pi_shift_atom_photon_state.png)
Using now the photonic basis state $\ket{\uparrow^p_x}=\frac{1}{\sqrt{2}}(\ket{\uparrow^p} + \ket{\downarrow^p})$ and $\ket{\downarrow^p_x}=\frac{1}{\sqrt{2}}(\ket{\uparrow^p} - \ket{\downarrow^p})$ this pi shift would be the same as using a `CNOT`-gate.

### Math
![math_for_cnot_photon_atom_gate.png](assets/imgs/math_for_cnot_photon_atom_gate.png)

#### Measurements

First a state tomogrophy was performed, where we have an atom either in the $\ket{F,m = 1,1}$ or $\ket{F,m=2,2}$ state, and our photon in the photonic basis states shown above for the shift to act as a CNOT gate. The results are the following:
![truth_table_single_atom_single_photon.png](assets/imgs/truth_table_single_atom_single_photon.png)

where we observe the right state with a probability of 86%.
Here we are predominantly limited by three effects:

1. Optical mode matching, because the transverse overlap between teh free-space mode of the photon and cavity mode is 92(3)%
2. Quality of the preparation of the state $\ket{\uparrow^a}$, which is successful with a 96(1)% probability
3. Stability of the difference between the cavity resonance and teh frequency of the impinging laser pulse, which is about 300 kHz

The decisive feature that distinguishes a quantum gate from a classical one is the generation of entangled states from seperable input states. Here we prepare the following state which should result in the following maximally entangled bell state after the gate:
$\ket{\downarrow^a_x \downarrow^p_x} \rightarrow \ket{\Phi^+_{ap}} = \frac{1}{\sqrt{2}}(\ket{\downarrow^a \downarrow^p_x} + \ket{\uparrow^a \uparrow^p_x})$
![bell_state_generation_quantum_gate.png](assets/imgs/bell_state_generation_quantum_gate.png)

Limiting us are mostly:

1. Above-mentioned spatial mode mismatch between cavity and impinging photon (8(3)%)
2. Quality of our atomic state preparation, rotation and readout (5(1)%)
3. Imperfections in the photonic polarization measurement and detector dark counts (3%)
4. Small probability of having more than one photon in the impinging laser pulses (2%)

<!-- TODO: Finish writing notes about this paper -->
