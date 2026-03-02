README

This repository contains the code for reproducing all simulation results, supplementary results, and real data analyses reported in the paper. The folder structure and instructions are described below.

1. Main Simulation Studies (Cases 1–3)

The folders "Case1", "Case2", and "Case3" correspond to the three simulation settings presented in the main text.

Within each case folder:

"/DW" implements the method from Dette and Wu (2026).

"/New" implements the proposed method.

To reproduce the simulation results for a given case, run the "main.R" file in each corresponding subfolder.

For example:
Run "/Case1/DW/main.R" and "/Case1/New/main.R" to reproduce the results reported in Table 1.

Similarly, running the corresponding "main.R" files under Case2 and Case3 reproduces the remaining simulation results in the main text.

2. Supplementary Simulation Studies (Cases 4–5)

The folders "Case4" and "Case5" correspond to the additional simulation settings reported in the Supplementary Materials.

To reproduce the results:

Run "/Case4/main.R" to obtain the results in Table S.1.

Run "/Case5/main.R" to obtain the results in Table S.2.

3. SCB for Marginal Mean Functions (Supplementary Section 4.3)

The folder "marginal_mean" contains the code for Supplementary Section 4.3: SCB for Marginal Mean Functions.

Subfolders include:

"/BCCK" — method from Belloni et al. (2015).

"/DW" — method from Dette and Wu (2026).

"/New" — the proposed method.

To reproduce the results in Tables S.3–S.4, run the "main.R" file in each subfolder.

4. Computational Efficiency (Supplementary Section 4.4)

The folder "time_consuming" contains the code for Supplementary Section 4.4: Computational Efficiency.

This section compares the computational cost of the proposed method with the method from Reference 1.

To reproduce Figure S.4, run: "/time_consuming/main.R"

5. Real Data Analyses

(1) Sydney Temperature

The folder "Sydney_temperature" contains the real data analysis of average temperature in Sydney.

Run "/Sydney_temperature/main.R" to reproduce all results reported in Section 7.1 of the main text.

(2) Australian Fertility Rate

The folder "Fertility" contains the real data analysis of Australian fertility rates.

Run "/Fertility/main.R" to reproduce all results reported in Section 7.2 of the main text.

(3) Montreal Temperature

The folder "Montreal_temperature" contains the additional real data analysis of temperature data from the Montreal region.

Run "/Montreal_temperature/main.R" to reproduce all results reported in Section 3 (Additional Real Data Analysis) of the Supplementary Materials.

Notes

Each "main.R" file is self-contained for its corresponding experiment.

If you encounter any issues, ensure that the mentioned R packages are installed. 

A random seed is specified in each script to reproduce the manuscript and supplemental results; 
slight variations may occur due to hardware differences, and the presented results have been reorganized from the computed outputs.

Simulation time depends on hardware; on a system with an 11th Gen Intel Core i7 CPU (2.80 GHz) and 16 GB RAM,
a full run takes around 24 hours or more. It's normal for the code to run without immediate feedback.
