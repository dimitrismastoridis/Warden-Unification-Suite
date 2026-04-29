# Warden Sector Precision Unification and Monte Carlo Convergence Suite
# U(4) Precision Gauge Coupling Unification
[![DOI](https://zenodo.org/badge/1213597114.svg)](https://doi.org/10.5281/zenodo.19895080)
This repository contains the computational framework and numerical proofs for precision gauge coupling unification within the $U(4)$ Grand Unified Theory. 

The codebase utilizes  3-loop Renormalization Group Equations (RGEs) to map the trajectories of the Standard Model gauge couplings ($\alpha_1, \alpha_2, \alpha_3$) from the electroweak scale up to the unification scale. It demonstrates that under the topological constraints of the $U(4)$ vacuum, the forces achieve flawless geometric unification and rigorously satisfy an exact integer quantization boundary condition.

## Overview of the Physics

Unlike standard $SU(5)$ or $SO(10)$ models that rely on arbitrary supersymmetric mass spectra to force unification, the $U(4)$ framework relies on deterministic topological transitions:
1. **The Warden Threshold (8.2 TeV):** The onset of new scalar degrees of freedom that flatten the running of the weak coupling, preventing it from overshooting.
2. **The Melt Scale (258 TeV):** The threshold of the Cho-Savvidy background. At this scale, a phenomenologically extracted Transverse Decoupling vector ($\Theta_D$) acts as a non-perturbative mass-gap friction on the couplings.
3. **The Spin-Transition Shift:** The unified couplings intersect in the $\overline{MS}$ scheme at $\approx 37.331$. A topological spin-shift ($+44/3\pi$) transports the unified couplings into the native topological vacuum, locking flawlessly onto **Target 42** ($\alpha_{native}^{-1} = 42.000$).

## Repository Structure

The framework is divided into two strict, uncoupled scripts to maintain absolute epistemological integrity. **Zero parameters are forced during execution.** ### 1. `monte_carlo_engine.py` (The Blind Predictor)
This script tests the physical robustness of the theory. It takes the frozen $\Theta_D$ array and subjects the Standard Model experimental inputs ($M_Z$, $\alpha_{em}$, $M_{top}$, $\alpha_s$) to random Gaussian noise.
* **Method:** Runs a 500-universe Monte Carlo simulation completely blind.
* **Result:** Achieves a 100% geometric convergence rate.
* **Prediction:** Organically predicts a low-energy weak mixing angle of $\sin^2\theta_W = 0.23129 \pm 0.00021$, perfectly matching current PDG experimental data.
================================================================
     U(4) RAW BLIND PREDICTIONS
================================================================
Convergence Rate     : 100.0%
OUTPUT sin^2 theta_W : 0.23129 ± 0.00021
OUTPUT M_GUT Scale   : 3.217e+16 GeV
================================================================

### 2. `unification_plotter.py` (The High-Precision Visualizer)
This script uses the exact predicted parameters from the Monte Carlo engine and applies ultra-strict mathematical tolerances (`rtol=1e-11`, `atol=1e-11`) to solve the stiff differential equations.
* **Method:** Maps the exact mean trajectories of the universe.
* **Result:** Demonstrates a flawless "zero-triangle" intersection of all three forces at $M_{GUT} \approx 3.217 \times 10^{16}$ GeV.
* **Output:** Generates the publication-ready figure `unificationplot.png`.
================================================================
MS_bar Intersection: 37.331
Native Target      : 42.000
Plot generated successfully as 'unificationplot.png'
================================================================
## Installation & Requirements

The scripts are written in Python 3. To run the RGE solvers and generate the plots, you will need the standard scientific computing stack.

```bash
pip install numpy scipy matplotlib

