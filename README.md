# Morphological Decomposition and Bar Identification of galaxies in the TNG100 simulation

**[ WORK IN PROGRESS: Analysis extension to z=6 underway, paper in prep. ]**

This repository hosts the codebase used for the Master's Thesis in Astronomy and Astrophysics at Sapienza University of Rome and the extensions for the paper in prep.

The project focuses on the kinematic analysis and morphological classification of galaxies extracted from the **IllustrisTNG (TNG100-1)** cosmological simulation. The code implements advanced decomposition techniques to disentangle galactic structural components and characterize non-axisymmetric features such as stellar bars across cosmic time.

## Scientific Methodology

The analysis framework is built upon two main theoretical pillars aimed at dissecting the internal structure of simulated galaxies:

### 1. Kinematic Decomposition (MORDOR)
The Python-based module adapts the MORDOR algorithm to the TNG100 volume. This method performs a dynamical decomposition of the galaxy based on the phase-space coordinates of stellar particles.

The core metric used is the **orbital circularity** ($\epsilon$), defined as:

$$\epsilon = \frac{j_z}{j_{\text{circ}}(E)}$$

where $j_z$ is the specific angular momentum component aligned with the galaxy's rotation axis, and $j_{\text{circ}}(E)$ is the angular momentum of a circular orbit at the same binding energy.

By analyzing the distribution of circularities, the code statistically identifies and separates distinct stellar components:
* **Disc components:** Associated with high circularity values ($\epsilon \sim 1$), representing rotationally supported structures.
* **Bulge/Spheroid components:** Associated with circularity values around zero ($\epsilon \sim 0$), representing pressure-supported structures with random motion.
* **Stellar Halo:** Identified based on energy and binding criteria.

This decomposition allows for the automated morphological classification of galaxies and the study of the transition from dispersion-dominated to rotation-dominated systems.

### 2. Bar Analysis and Fourier Decomposition
The C-based module (`Hfourier.c`) is designed to investigate the properties of non-axisymmetric perturbations within the identified stellar discs.

This tool computes the Fourier modes of the surface density distribution of the stellar particles. By calculating the amplitude of the Fourier components (specifically the $m=2$ mode), the code quantifies the strength of the non-axisymmetric features. This theoretical approach is used to:
* Confirm the presence of stellar bars.
* Measure the strength and length of the bar.
* Analyze the evolution of the bar structure in relation to the galaxy's dynamical history.

## Project Scope and Evolution

The current iteration of the code analyzes galaxy samples from the TNG100 simulation. Work is currently in progress to extend the analysis pipeline to high redshifts (up to $z=6$). This involves processing approximately 100 simulation snapshots to reconstruct the full evolutionary tracks of the progenitors and study the onset of morphological features in the early Universe.
