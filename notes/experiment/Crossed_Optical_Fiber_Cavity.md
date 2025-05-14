---
id: "Crossed_Optical_Fiber_Cavity"
date: "2025-05-14"
subtitle: 
tags:
    - cavity
    - fibers
    - experiment
---

# Crossed_Optical_Fiber_Cavity

At the heart of the `CavityX` setup lies the crossed optical fiber cavity (hence the name).
It consists of two fiber cavities arranged at and angle of 90 degrees to each other:

![crossed_fiber_cavity_setup_closeup_image.png](assets/imgs/crossed_fiber_cavity_setup_closeup_image.png)

The length of the two resonators is the inspiration behind their regularly used names: short (SC; red) and long cavity (LC; blue).
The coupling of light into and out of each resonator occurs via the same single-mode fiber, (a) and (b) since the cavities are designed to be single-sided.
This is achieved by coating the single-mode fiber surfaces with a low reflective dielectric mirror coating and the multi-mode fiber surfaces, (c) and (d), with a highly reflective dielectric mirror coating.

A requirement is that the in- and outcoupling channel of the optical cavities requires the transverse mode of the optical field to be preserved. This is why single-mode fibers with a highly transmissive mirror coating are used, whereas on the other hand you have the opposing multi-mode fibers which are machined with a highly reflective coating.
Using multi-mode fibers is advantageous because:
    - Multi-mode fibers have large core diamater -> allows for high and robus light coupling efficiency -> latter is beneficial for cavity transmission measurements
    - Produce better concave surface structures when a $CO_2$ machining process is applied

## Cavity Parameters

<!--TODO: get parameters from Gianvito -->

---

`Note: Short cavity cannot be coupled to the stronged $D_2$-line transition $F=2,m_F=\pm2 \leftrightarrow F'=3,m_F=\pm3$ since it only supports linearly polarized cavity modes.`

---

Furthermore, we observe the unequal radii of curvatures (ROCs), which has been chosen intentionally in order to move the cavity mode waist towards the in-/outcoupling mirror, thus increasing fiber-cavity mode matching.

## Mounting Of The Fiber Cavities

The goal here is to mount very fragile fibers securely and machanically decoupled from the environment, while at the same time allowing them to be optically accessible and be able to move htem over a broad range with sub-micrometer precision.

One way to achieve this is to use UHV compatible glue in order to mount the fibers onto piezo elements directly, which unfortunately bears the risk of the fiber moving unpredictably during the curing of the glue.

The chosen method of mounting the fibers is shown in the following picture:
![cavityX_fiber_mount.png](assets/imgs/cavityX_fiber_mount.png)

`a:` An early version of the fiber cavity holder. The concept of this mount has been adapted to the present version. 
However, the dimensions of the constituting elements were adjusted to fit into the interior of the vacuum chamber. 
(1), (3) and (4) are custom-designed Macor elements. 
(2) is a shear piezo soldered to two copper wires. 
Titanium screws and Inconel springs (5) fix the upper Macor plate (4) to the lower one (3).
The optical fiber (6) is clamped by the Macor elements (3) and (4).
`b:` The fiber cavity holder installed into the vacuum chamber. 
The labels of subﬁgure (a) are adopted. One piezo slip-stick translation stage (7) can be observed: it enables the vertical movement of the optical fiber. 
(8) labels an in-vacuum mirror that reflects a magneto-optical trap beam. 
The dashed rectangle surrounds the position of the optical fiber from which on the copper coating has been etched away. 
It is necessary to have a coating-free optical ﬁber at the end facet for the CO2 beam machining process of the concave mirror surfaces.


##### References
Niemietz, D. Nondestructive detection of photonic qubits with single atoms in crossed ﬁber cavities. (Technische Universität München,2021), https://mediatum.ub.tum.de/1601857
