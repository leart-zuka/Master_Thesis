---
id: "fibres"
date: "2025-05-13"
subtitle: 
tags:
    - fibers
---

# Fibers

This is going to be a summary of informations regarding fibers

## Design

Optical fibers consist of two layers:
    - The core layer ($n_1$)
    - The cladding layer ($n_2$)

What differentiates these two layers is the difference in refractive indices between two ($n_2 < n_1$). Typically they're chosen so that there is total internal reflection between the two layers for minimal losses.

## Types of Fibers

There are 3 main different types of fibers:
    - Single-Mode Fibers
    - Polarization-Maintaining Fibers (PM Fibers) ![pm_fiber.png](assets/imgs/pm_fiber.png)
    - Multi-Mode Fibers 

What differentiates a multi-mode fiber and a single-mode fiber is the diameter of the core layer, as a bigger core allows for modes to enter the fiber and form these Hermite-Gaussian modes:
![difference_single_mode_multi_mode.png](assets/imgs/difference_single_mode_multi_mode.png)

### Notes

The reason why one would never use PM fibers over single-mode fibers in Quantum-Teleportation protocol experiments, is because you need that extra degree of freedom in the form of polarization in these experiments.
Using a PM fiber would destroy your qubit state, thus ruining the experiment
