---
id: "quantum_teleportation_with_one_mediating_photon"
date: "2025-04-29"
subtitle: 
tags:
    - paper
---

# quantum_teleportation_with_one_mediating_photon

This is going to be my summarization of the paper `Quantum Teleportation between Remote Qubit Memories with Only a Single Photon as a Resource` 
(DOI:[10.1103/PhysRevLett.126.130502](https://doi.org/10.1103/PhysRevLett.126.130502))

## Intro

This paper showcases a new protocol that performs unconditional teleportation without the use of a preshared entangled qubit pair, but rather through the use of only a single photon as an `ex-ante`(latin: "before the event") prepared source. Because what you usually do is you have a qubit on one side, which you then entangle with e.g. a photon. And now you want to copy the state of the first qubit to the another qubit "on the other side". In order to do so you also entangle the target qubit with e.g. a photon, send both photons (from source and target) through a lossy channel (e.g. a fibre) and perform a bell state measurement (BSM). The result from the bell state measurement then gives you the information on what feedback you have to send to the target qubit in order to reproduce the state from the source qubit. That feedback comes in the form of e.g. X/Y/Z gates. A classical circuit looks like this:
![classical_quantum_teleportation.png](assets/imgs/classical_quantum_teleportation.png)
However there are two main challenges residing in this design:
    1. Source must keep qubit alive while entanglement is generated over the quantum channle
    2. Entanglement distribution over the lossy channel must be heralded (in other words: Alice should only perform bell state measurement and all the stuff if a photon was actually detected)

## New design

This new design now in principle, allows for unconditional teleportation without the necessity of the commonly employed pre-shared entanglement resource, since we're only using one photon travelling from Bob to Alice.
What's also nice is, if the photon is lossed over the course of travelling through the lossy channel to Alice, Alice's qubit isn't affected. This allows us to repear the protocol until a successful photon transmission is heralded.
What's also interesting is that the entanglement is generated on the fly and is in principle deterministic. Furthermore in this scheme a BSM is implemented, with no fundamental efficiency limitation. 
Detecting the photon at Alice's side let's photon act as two way herald:
1. Entanglement generation on Bob's side (since before it went to Alice it bounced off of Bob's node, creating entanglement; detecting the photon on Alice's side means that we still have that entangled photon)
2. BSM on Alice's side (since it arrived on Alice's side, we can use it for BSM)

This protocol practically achieves high teleportation reates and fidelities, and is ideally suited for future quantum networks
This is what it looks like:
![new_quantum_teleportation_design_only_one_photon.png](assets/imgs/new_quantum_teleportation_design_only_one_photon.png)

## Protocol steps
1. Trap atoms @ both nodes
2. Send light to NPBS @ target loc. to reflect off of target node (acts as CNOT between target qubit & photon => light changes from $\ket{A}$ to $\ket{D}$)
3. Reflect photon off of source node (CNOT between source qubit & photon => change polarization again [quantum_gate_between_flying_photon_and_trapped_atom](papers/quantum_gate_between_flying_photon_and_trapped_atom.md)) 
4. Perform $\frac{\pi}{2}$ rotation around y-axis to Alice & measure in normal basis (here z-basis)
5. Measure photon in polarization basis
6. Send classical info back to Bob

| Heralded State (photon + Alice's atom) | Feedback to Bob     |
|----------------------------------------|----------------------|
| \|A⟩                                   | none                 |
| \|D⟩                                   | \(R_{x}(\pi)\)       |
| \|↑⟩                                   | Z-Gate               |
| \|↓⟩                                   | none                 |

### Math
![bsp_for_new_teleporation_scheme.png](assets/imgs/bsp_for_new_teleporation_scheme.png)
![bsm_and_feedback_for_new_teleporation_scheme.png](assets/imgs/bsm_and_feedback_for_new_teleporation_scheme.png)

#### Measurements

Two $^{87}Rb$ atom trapped in two high-finesse optical cavities connected via a 60m long single mode optical fiber. Each atom carries one qubit information
encoded in the states $\ket{5^2 S_{\frac{1}{2}}, F=2, m_F=2}:=\ket{\uparrow_z}$ and $\ket{5^2 S_{\frac{1}{2}}, F=1, m_F=1}:=\ket{\downarrow_z}$. Both cavities
are actively tuned to the atomic resonance $\ket{\uparrow_z}\leftrightarrow \ket{e}:=\ket{5^2 P_{\frac{3}{2}}, F'=3, m_F=3}$.
We use a weak coherent pulse (average photon number $\langle n \rangle=0.07$.
Excluding optical pumping ($200 \mu s$), the entire protocol takes $25.5 \mu s$ (duration of Raman pulse: $4\mu s$ for a pi/2 pulse).
Repitition rate of the experiment is 1kHz, probability for successful transmission is 8.4% giving us a teleportation rate of 84Hz (taking the weak coherent pulses into account we get a total rate of 6 Hz).
In order to check the performance of the protocol prepare Alice's qubit in six mutually unbiased states and teleport these to Bob, and then measure Bob's state in order to get a quantum state tomogrophy:
![one_photon_teleportation_result.png](assets/imgs/one_photon_teleportation_result.png)

We are limited by namely:
1. Qubit decoherence (6.0%)
2. 2-photon contributions in the laser pulses (3.9%)
3. State-preparation (1.4%)

Scanning over the mean photon number gives us the following result:
![one_photon_teleportation_mean_photon_number_result.png](assets/imgs/one_photon_teleportation_mean_photon_number_result.png)
Here we get the highest values for a vanishing mean photon number and a decrease due to higher photon-number contributions. Classical threshold is at a photon number of around 1.0, giving us a fidelity of 2/3.

Increasing delays between the gates and classical feedback to simulate a larger distance between two networks gives us following circuit:
![one_photon_teleportation_simulated_delay_circuit.png](assets/imgs/one_photon_teleportation_simulated_delay_circuit.png)
The reason for the placement of the delays is:

- Between state-preparation pulses to simulate a communication after Alice's qubit is prepared. 
- Afterwards another longer one to mimic propagation time of the photon in the fiber
- A third one to take into account the propagation time of the two feedback signals.

This gives us the following results
![one_photon_teleportation_simulated_delay_result.png](assets/imgs/one_photon_teleportation_simulated_delay_result.png)
A delay of $40 \mu s$ simulates a distance of 8km. It should be noted however that this is only a simple simulation and doesn't account for the problems with actually using fibers this long.  
