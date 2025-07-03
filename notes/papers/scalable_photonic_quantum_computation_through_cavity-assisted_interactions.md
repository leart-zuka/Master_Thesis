---
id: "scalable_photonic_quantum_computation_through_cavity-assisted_interactions"
date: "2025-05-27"
subtitle: 
tags:
    - paper
---

# scalable_photonic_quantum_computation_through_cavity-assisted_interactions

This is going to be my summarization of the paper `Scalable Photonic Quantum Computation through Cavity-Assisted Interactions` 
(DOI: [10.1103/PhysRevLett.92.127902](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.92.127902))

## Intro

The paper introduces a new scalable quantum computation scheme, which encodes informations in the polarization of photonic qubits, through interactions with a trapped atom inside of an optical cavity.

This scheme can then be used to achieve controlled gate operations and also act as a high efficiency single-photon detector.

It can then be scaled up to many qubits, and be integrated with protocols for the realization of quantum networks.


`Advantage of using photonic qubits for quantum computation:

- Can be readily scaled up by generating many single-photon pulses

Disadvantage of using photonics qubits for quantum computation:

- Quantum gates between single-photon pulses very hard :( (photon-photon coupling rate too small to allow for meaningful operations)`

Therefor, use single-photon detector (aka an atom) to achieve non-linear interaction between photons through feed forward. 

**IMPORTANT**: Required efficiency $\alpha$ needs to be very high in order to be viable (gate success probability $p\simeq0.99, \alpha\geq0.999987$)

In this scheme: Combine photonic qubits with power of strong atom-photon coupling through cavity. Cavity in far-off resonant interactions provides effective Kerr nonlinearity for the input light.

New protocol has following advantages:

1. Different interaction mechanism has much larger effective interaction rate, allowing for implementation of CNOT gate
2. Phase flip is insensitive to variation of atom-photon coupling rate
3. Pulse shapes for pairs of interacting single photons suffer small changes
4. Noise properties of the scheme are very favorable

## Actual Scheme

Photonic qubit is encodid in the polarization states: $\ket{h}$ and $\ket{v}$

Single qubit operations are easy, because we can do polarization rotations, aka just use waveplates.
Problems are the nontrivial-two-qubit interactions. Here we're going to look at the controlled phase-flip (CPF) gate, which `for the two qubits i and j flips the phase of the input if both qubits are in` $\ket{h}$. 
Using the CPF gate plus normal single qubit rotations, realize universal quantum computation.

### Setup

![cpf_gate_scheme.png](assets/imgs/cpf_gate_scheme.png)

Implementation of CPF gate between two photons is implemented by "bouncing" them off of the cavity after one another, with an atom being trapped inside of the cavity.

Atom has three relevant levels, and is initialized in the superposition of two ground states $\ket{\Phi_{ai}} = \frac{\ket{0}+\ket{1}}{\sqrt{2}}$.

Atomic transition $\ket{1} \rightarrow \ket{e}$ is resonant to the h polarization of our photon, and is also coupled to the cavity mode $a_h$.

There are two critical steps which make up the scheme:

1. CPF gate between the photon and the atom is created by just bouncing a photon j off of the cavity and is described by: $U_{aj}^{CPF} = e^{i\pi\ket{0}_a\bra{0}\otimes\ket{h}_j\bra{h}}$


The CPF gate looks like this btw: ![cpf_gate_matrix_notation.png](assets/imgs/cpf_gate_matrix_notation.png)

2. Composition of CPF gates between atom and pulses j and k, describes unitary operator: $U_{jk}^{CPF} = e^{i\pi\ket{h}_j\bra{h}\otimes\ket{h}_k\bra{h}}$, while restoring atom back to $\ket{\Phi_{ai}}$

When photon is *v* polarized, it will just get deflected and not experience a phase shift.

Is the photon *h* polarized, it is resonant with the bare cavity mode if atom is in $\ket{0}$. It will then aquire a phase shift of $e^{i\pi}$.

If the atom is in state $\ket{1}$ however, frequency of dressed cavity mode is significantly detuned from freq. of incoming pulse, thus pulse will get reflected as if it was reflected off of a mirror.

Following steps then make up our CPF gate between our two pulses:

1. Reflect photon j from cavity -> $U_{aj}^{CPF}$
2. Apply $\pi/2$ rotation on atom
3. Reflect photon k from cavity -> $U_{ak}^{CPF}$
4. Apply $\pi/2$ rotation on atom
5. Reflect photon j from cavity again -> $U_{aj}^{CPF}$

giving us following unitary operation:

$U_{jk}^{CPF} \ket{\psi_{jk}}\otimes\ket{\Phi_{ai}} = U_{aj}^{CPF}\cdot R_a(\pi/2)\cdot U_{ak}^{CPF}\cdot R_a(\pi/2)\cdot U_{aj}^{CPF} \ket{\psi_{jk}}\otimes\ket{\Phi_{ai}}$

#### QND 

Can also use this scheme to implement a `quantum nondemolition (QND)` measurement of the photon number in the pulse.

For this we need to:

1. Prepare atom in $\ket{\Phi_{ai}}$
2. Reflect photon pulse from cavity
3. Apply a $R_a(\pi/2)$
4. Measure atomic state in basis $\{\ket{0},\ket{1}\}$

If outcome is "0" if pulse has a photon in $\ket{h}$.

## Theoretical Model

Initial state of a pulse can be expressed as $\ket{\psi_p}_j = c_{hj}\ket{h}_j + c_{vj}\ket{v}_j$.

$c_{hj}$ and $c_{vj}$ are superposition coefficients.

The polarization eigenstates $\ket{h}_j$ and $\ket{v}_j$ have the form $\ket{h;v}_j = \int_0^T f_j(t)a_{h;v}^{in\dagger}(t)dt \ket{vac}$.

$f_j(t)$ is the `normalized pulse shape` as a function of time t, T is the `pulse length duration`, and $a_{h;v}^{in}$ are the `one-dim. field operators (cavity input operators)` with standard commutator relationship. 
$\ket{vac}$ is just the vacuum of all optical modes.

Cavity mode $a_h$ is driven by corresponding cavity input operator $a_h^{in}(t)$ via:
![cavity_mode_driving.png](assets/imgs/cavity_mode_driving.png)

$\kappa$ is the cavity decay rate and hamiltonian is:

$H = \hbar g (\ket{e}\bra{1}a_h + \ket{1}\bra{e}a_h^\dagger)$ (3)

Hamiltonian describes coherent interaction between atom and cavity mode $a_h$.

Detuning $\laplace$ in (2), is 0 for our scheme. 

Cavity output $a_h^{out}(t)$ is connected to input via `input-output` relation:

$a_h^{out}(t) = a_h^{in}(t) + \sqrt{\kappa}a_h$ (4).

---
NOTE: For *v* component of the pulse, we simply have $a_v^{out}(t) = a_v^{in}(t)$ as it just gets reflected.

---

Equations (2)-(4) determine evolution of joint state of atom and photon.

Analytically we can start off by saying, if we start off with out atom in $\ket{0}$, Hamiltonian does not play a role, ergo we find

$a_h^{out}(t) \approx \frac{i\laplace - \kappa/2}{i\laplace + \kappa/2} a_h^{in}(t)$.

Here high freq. components of our field operators have been discarded (valid if pulse shapre changes slowly over time compared to cavity decay rate).

Under this approx. we have, $a_h^{out}(t) \approx - a_h^{in}(t)$ for resonant interaction ($\laplace = 0$), so *h* component acquires $\pi$ phase.

If atom is in $\ket{1}$ however, we need to take hamiltonian into account, where then for **strong coupling**, we need to remember that our two dressed cavity modes, have frequencies that are actively detuned from input pulse by$\laplace = \pm g$. For case $g \gg \kappa$ we get $a_h^{out}(t) \approx a_h^{in}(t)$, ergo no phase shift.

This confirms that we actually have a CPF gate.

## Simulation Results

Atomic spont. emission noise is described by an imaginary part
$(-i\gamma_s/2)(\ket{e}\bra{e} - \ket{1}\bra{1})$ 
in the Hamiltonian H, with $\gamma_s$ being the spontaneous emission rate from state $\ket{e}$.

Input pulse shape is assumed to be gaussian with the following shape $f_j(t) \propto exp(-(t-T/2)^2/(T/5)^2)$, where t goes from 0 to T.

### Input Output Shape

![cpf_gate_simulation_results_paper.png](assets/imgs/cpf_gate_simulation_results_paper.png)

2a shows the shape function for the input pulse (solid curve) and the reflected pulse with the atom being in $\ket{0}$ (dashed curve) and $\ket{1}$ (dotted curve). The simulation shows that the output pulse basically has teh same shape as the input pulse if the pulse duration is $T \gg 1/\kappa$.

2b shows the gate fidelity F for different pulse durations T in units of $\kappa$.
For $T = 240/\kappa$, corresponding to a pulse width of $T/5 \sim 1\mu s$, the gate fidelity is $99.9\%$. Furthermore the shape of the output pulse is also very insnsitive to variation of the coupling rate g in the typical parameter region (change is smaller than $10^{-4}$ for g varying from $\kappa$ to $6\kappa$).

Dominant noise arises from atomic spontaneous emission, leading to vacuum-state output. This noise is a `leakage error` which means it is outside of the qubit Hilbert space.

2c shows the Probability $P_s$ of spont. emission loss as a function of $g/\kappa$ for the input state $\ket{1} \otimes \ket{h}$, assuming $\gamma_s=\kappa$. Curve is well simulated by $P_s \approx 1/(1+2g^2/\kappa\gamma_s)$.

##### References
L.-M. Duan and H.J. Kimble - Scalable Photonic Quantum Computation through Cavity-Assisted Interactions (2004), https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.92.127902
