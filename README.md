# Warden Sector Precision Unification and Monte Carlo Convergence Suite


[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview
This repository contains the complete numerical framework used to verify the exact gauge coupling unification presented in the paper: **"Topological Vacuum Misalignment and the 8.2 TeV Warden Resonance: A Path to Gauge Coupling Unification"** by D. Mastoridis, K. Kalogirou, and P. Razis.

The code evaluates the Renormalization Group Equations (RGEs) up to the three-loop level, incorporating the topological melting of the Warden condensate at 258 TeV and the subsequent activation of four complex scalar bi-doublets. It includes the deterministic $\overline{MS}$ phase transition mechanics and a high-precision Monte Carlo suite to propagate experimental uncertainties.

## Physics Highlights
* **High-Energy RGEs:** Simulates 1-loop, 2-loop, and 3-loop gauge interactions, alongside top-quark $\overline{MS}$ Yukawa corrections.
* **Strict Numerical Stability:** Utilizes implicit Runge-Kutta (`Radau`) integration with microscopic tolerances for stiff ODEs. The strict implementation includes dynamic event-termination to rigorously reject non-physical trajectories (Landau poles and Yukawa blowups).
* **Topological Thresholds:** Applies exact geometric ($\Theta_S$) and transverse decoupling ($\Theta_D$) matching conditions at the 258 TeV Sphaleron barrier.
* **Monte Carlo Verification:** Performs bounded simulations drawing from Particle Data Group (PDG) distributions to verify the 100% concurrent unification rate at $M_{GUT} \approx 3.23 \times 10^{16}$ GeV.
* **Precision Post-Diction:** Confirms the derivation of the weak mixing angle ($\sin^2\theta_W \approx 0.23129$).

## Repository Structure
This suite includes three main computational scripts:
* `monte_carlo_standard.py`: The baseline statistical suite utilizing standard Gaussian error propagation for the $U(4)$ framework.
* `monte_carlo_strict.py`: The advanced statistical suite. Executes rigorous uncertainty propagation using $3\sigma$ truncated normal distributions and strict RGE event-termination to ensure absolute physical validity across all simulated universes.
* `plot_unification.py`: Generates the high-resolution, publication-ready visualization of the exact gauge coupling intersection , applying the final geometric scaling to reach $\alpha_{native}^{-1} \equiv 41.9$.

## Dependencies
The suite relies on standard scientific Python libraries for integration, root-finding, and plotting. To install:
```bash
pip install numpy scipy matplotlib

Usage and Expected Outputs
1. The Standard Monte Carlo Suite
Run the baseline 5,000-universe statistical analyzer to evaluate the parameter space using standard Gaussian error propagation:

Bash
python monte_carlo_standard.py
Expected Output:

Plaintext
Initiating 3-Loop Monte Carlo Predictor (5000 Universes)...

================================================================
     U(4) 3-LOOP PREDICTIONS WITH FULL ERROR PROPAGATION
================================================================
Convergence Rate     : 100.0% (5000/5000 runs unified)
OUTPUT sin^2 theta_W : 0.23129 ± 0.00021
OUTPUT M_GUT Scale   : 3.228e+16 ± 1.668e+15 GeV
================================================================
2. The Strict Event-Terminated Suite
Run the highly rigorous analyzer. This version uses truncated normal distributions and actively tracks the RGE flow, immediately terminating any universe that develops a Landau pole or Yukawa blowup:

Bash
python monte_carlo_strict.py
Expected Output:

Plaintext
Initiating Strict 3-Loop Monte Carlo Predictor (500 Universes)...
Running... 500/500 complete.

================================================================
     U(4) STRICT 3-LOOP PREDICTIONS WITH EVENT TERMINATION
================================================================
Convergence Rate     : 100.0% (500/500 runs unified)
OUTPUT sin^2 theta_W : 0.23128 ± 0.00020
OUTPUT M_GUT Scale   : 3.238e+16 ± 1.594e+15 GeV
================================================================
3. Generate the Unification Plot
To generate the 3-loop unification trajectory graph using the central PDG values:

Bash
python plot_unification.py
This will render the graph and save it to your directory as FigureW_GUT_Unification.pdf.

## Dependencies
The suite relies on standard scientific Python libraries. To install:
```bash
pip install numpy scipy
