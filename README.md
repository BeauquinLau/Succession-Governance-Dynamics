# Tipping Points and Resilience: An Evolutionary Game Theory of Governance Choice in Family Firm Succession

This repository contains the MATLAB source code for the numerical experiments presented in the associated manuscript. This document provides detailed instructions to ensure the clear and accurate reproduction of all figures.

## Repository Structure

The repository is organized into two main subdirectories. The code is structured to mirror the analyses in the paper.

```
Simulation/
│
├── Figures_pdf/
│   ├── S1.pdf, S1_d.pdf, S2.pdf, S2_d.pdf
│   ├── Dynamics_n_1.pdf, Dynamics_n_2.pdf
│   ├── Dynamics_phi_NPM_1.pdf, Dynamics_phi_NPM_2.pdf
│   ├── Dynamics_f_gamma_1.pdf, Dynamics_f_gamma_2.pdf
│   ├── A.pdf, B.pdf, C.pdf, D.pdf, E.pdf, F.pdf
│   ├── G.pdf, H.pdf, I.pdf, J.pdf, K.pdf, L.pdf
│   └── M.pdf, N.pdf, O.pdf, P.pdf, Q.pdf, R.pdf
│
└── MATLAB_code/
    │
    ├── System_Dynamics_Validation/
    │   ├── S1.m, S1_d.m
    │   └── S2.m, S2_d.m
    │
    └── Sensitivity_Analysis/
        │
        ├── Analysis_n/
        │   ├── Dynamics_n_1.m, Dynamics_n_2.m
        │   └── Trajectories/
        │       └── A.m, B.m, C.m, D.m, E.m, F.m
        │
        ├── Analysis_phi_NPM/
        │   ├── Dynamics_phi_NPM_1.m, Dynamics_phi_NPM_2.m
        │   └── Trajectories/
        │       └── G.m, H.m, I.m, J.m, K.m, L.m
        │
        └── Analysis_f_gamma/
            ├── Dynamics_f_gamma_1.m, Dynamics_f_gamma_2.m
            └── Trajectories/
                └── M.m, N.m, O.m, P.m, Q.m, R.m
```

## System Requirements

-   **MATLAB**: The code has been tested on MATLAB R2025a and is expected to be compatible with recent versions.
-   **Toolboxes**: No special toolboxes are required. The code relies only on the base MATLAB environment.

## Instructions for Reproduction

This repository provides executable scripts to generate all figures from the manuscript's numerical exploration section.

> **General Instruction:** To run any script, please first set MATLAB's current working directory to the script's location. All generated figures will be saved as `.pdf` files in the same directory.

### 1. Dynamic Evolution of Governance Configurations (Section 5.1)

These scripts validate the main theorems by illustrating the system's phase portraits.

-   **Location**: `MATLAB_code/System_Dynamics_Validation/`
-   **`S1.m`**: Generates the trajectory visualization for the baseline cooperative model.
    -   **Output**: `S1.pdf` (Corresponds to Figure 1(a)).
-   **`S1_d.m`**: Generates the regional location diagram for the baseline model.
    -   **Output**: `S1_d.pdf` (Corresponds to Figure 1(b)).
-   **`S2.m`**: Generates the trajectory visualization for the extended conflictual model.
    -   **Output**: `S2.pdf` (Corresponds to Figure 2(a)).
-   **`S2_d.m`**: Generates the regional location diagram for the extended model.
    -   **Output**: `S2_d.pdf` (Corresponds to Figure 2(b)).

📌 **To run (example for Figure 1):**

```matlab
% In MATLAB, navigate to the 'System_Dynamics_Validation' directory
run('S1.m');
run('S1_d.m');
```

### 2. Sensitivity Analysis of Exogenous Parameters (Section 5.2)

These scripts demonstrate how the system's tipping point (saddle point $E_5$) and evolutionary paths shift as key parameters are varied.

#### 2.1. Impact of Monitoring Intensity (`n`)

-   **Location**: `MATLAB_code/Sensitivity_Analysis/Analysis_n/`

-   **Saddle Point Drift (Figure 3):**
    -   **Script**: `Dynamics_n.m` generates the diagrams showing the drift of the saddle point $E_5$ as `n` varies.
    -   **Outputs**: `Dynamics_n_1.pdf` (Fig. 3(a)) and `Dynamics_n_2.pdf` (Fig. 3(b)).

-   **Convergence Trajectories (Figure 4):**
    -   **Location**: `Trajectories/` subfolder.
    -   **Scripts**: `A.m` to `F.m`.
    -   **Outputs**: `A.pdf` (Fig. 4(a)) to `F.pdf` (Fig. 4(f)).

#### 2.2. Role of Agency Costs (`phi_NPM`)

-   **Location**: `MATLAB_code/Sensitivity_Analysis/Analysis_phi_NPM/`

-   **Saddle Point Drift (Figure 5):**
    -   **Script**: `Dynamics_phi_NPM.m` generates diagrams showing the drift of $E_5$ as `phi_NPM` varies.
    -   **Outputs**: `Dynamics_phi_NPM_1.pdf` (Fig. 5(a)) and `Dynamics_phi_NPM_2.pdf` (Fig. 5(b)).

-   **Convergence Trajectories (Figure 6):**
    -   **Location**: `Trajectories/` subfolder.
    -   **Scripts**: `G.m` to `L.m`.
    -   **Outputs**: `G.pdf` (Fig. 6(a)) to `L.pdf` (Fig. 6(f)).

#### 2.3. Influence of Affective Benefits (`f_gamma`)

-   **Location**: `MATLAB_code/Sensitivity_Analysis/Analysis_f_gamma/`

-   **Saddle Point Drift (Figure 7):**
    -   **Script**: `Dynamics_f_gamma.m` generates diagrams showing the drift of $E_5$ as `f_gamma` varies.
    -   **Outputs**: `Dynamics_f_gamma_1.pdf` (Fig. 7(a)) and `Dynamics_f_gamma_2.pdf` (Fig. 7(b)).

-   **Convergence Trajectories (Figure 8):**
    -   **Location**: `Trajectories/` subfolder.
    -   **Scripts**: `M.m` to `R.m`.
    -   **Outputs**: `M.pdf` (Fig. 8(a)) to `R.pdf` (Fig. 8(f)).

📌 **To run (example for the full analysis of parameter `n`):**

```matlab
% In MATLAB, navigate to the 'Sensitivity_Analysis/Analysis_n' directory
run('Dynamics_n_1.m');
run('Dynamics_n_2.m');

% Then navigate to the 'Trajectories' subfolder
cd('Trajectories');
run('A.m');
run('B.m');
% ... and so on.
```
