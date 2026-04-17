# Warden Sector Precision Unification and Monte Carlo Convergence Suite

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview
This repository contains the complete numerical framework used to verify the exact gauge coupling unification presented in the paper: **"Topological Vacuum Misalignment and the 8.2 TeV Warden Resonance: A Path to Gauge Coupling Unification"** by D. Mastoridis, K. Kalogirou, and P. Razis.

The code evaluates the Renormalization Group Equations (RGEs) up to the three-loop level, incorporating the topological melting of the Warden condensate at 258 TeV and the subsequent activation of four complex scalar bi-doublets. It includes the deterministic $\overline{MS}$ phase transition mechanics and a high-precision Monte Carlo suite to propagate experimental uncertainties.

## Physics Highlights
* **High-Energy RGEs:** Simulates 1-loop, 2-loop, and 3-loop gauge interactions, alongside top-quark MS-bar Yukawa corrections.
* **Numerical Stability:** Utilizes implicit Runge-Kutta (`Radau`) integration for stiff ODEs and `numpy.einsum` for exact 3-loop tensor contractions.
* **Topological Thresholds:** Applies exact geometric ($\Theta_S$) and transverse decoupling ($\Theta_D$) matching conditions at the 258 TeV Sphaleron barrier.
* **Monte Carlo Verification:** Performs 5,000 simulations drawing from Particle Data Group (PDG) Gaussian distributions to verify the 100% concurrent unification rate at $M_{GUT} \approx 3.228 \times 10^{16}$ GeV.
* **Precision Post-Diction:** Confirms the derivation of the weak mixing angle ($\sin^2\theta_W \approx 0.23129$).

## Dependencies
The suite relies on standard scientific Python libraries. To install:
```bash
pip install numpy scipy
