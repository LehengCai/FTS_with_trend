README

This repository contains the simulation and real data analysis code accompanying the paper "Inference for functional time series with time-varying trend". The folder structure and running instructions are described below. All paths are relative to the root folder "FTSSCR_code". Table, figure, and section numbers refer to the revised manuscript and its Supplementary Material.

1. Main Simulation Studies (Cases 1-3)

The folders "Case1", "Case2", and "Case3" correspond to the three simulation settings in Section 6 of the main text.

Within each case folder:

"DW" implements the method of Dette and Wu (2026).

"New" implements the proposed method.

For example, run "Case1/New/main.R" and "Case1/DW/main_new.R" for Case 1.

The output files contain empirical coverage probabilities and average SCR volumes over both the full domain and the interior domain. The results correspond to:

Case 1: Table 1 (full domain) and Table S.1 (interior domain).
Case 2: Tables S.2-S.3.
Case 3: Tables S.4-S.5.

The default settings are T = 400, 800, 1200, N = 100, and 500 Monte Carlo replications. The significance levels are alpha = 0.01, 0.05, 0.10, 0.20, corresponding to nominal coverage levels 1 - alpha.

Results are saved in the "result" subfolder as "T=...N=....csv" for the proposed method and "Wu_T=...N=....csv" for DW.
 


2. Supplementary Simulation Studies (Cases 4-5)

The folders "Case4" and "Case5" contain the additional simulations in Sections S.4.2-S.4.3 of the Supplementary Material.

Run "Case4/main.R" for the slowly decaying eigenvalue experiment associated with Table S.6. The data-generating process retains 100 functional principal components. The default sample sizes are T = 400, 800, 1200 and N = 100.

Run "Case5/main.R" for the small-sample coverage experiment associated with Table S.7. The default setting is (T, N) = (50, 400). The output reports the original coverage probabilities and those obtained after multiplying the SCR half-width by (1 + 20/T). 


3. Sensitivity Analyses (Supplementary Section S.4.4)

The folders "J_T", "J_N", "H", and "B" contain sensitivity analyses for the temporal knot number, within-curve knot number, block size, and number of Gaussian replications, respectively.

Run the following scripts for the corresponding experiments:

"J_T/New/main.R": Tables S.8-S.9; c0 = 0.6, 0.8, 1.0.

"J_N/New/main.R": Tables S.10-S.11; c0 = 1.2, 1.5, 1.8.

"H/New/main.R": Tables S.12-S.13; block-size multipliers 0.4, 0.5, 0.6.

"B/New/main.R": Tables S.14-S.15; Gaussian replication counts 500, 1000, 2000.

These scripts run one tuning value at a time, rather than automatically looping over all values. Repeat the corresponding run after changing the tuning value. Each output contains full-domain coverage probabilities and average SCR volumes and is saved in the local "result" subfolder.

For "J_T/New/main.R", change "c0" in the driver. 

For "J_N/New/main.R", change "c0" in the driver.

For "H/New/main.R", change the block-size multiplier in the "best_BB_index" assignment in "H/New/SCR.R" and update the matching "H=" filename label in the driver. The supplied multiplier is 0.6.

For "B/New/main.R", change the first argument of "mvrnorm()" in "B/New/SCR.R" and update the matching "B=" filename label in the driver. The supplied count is 2000.

4. SCBs for Marginal Mean Functions (Supplementary Section S.4.5)

The folder "marginal_mean" contains the experiments for the univariate mean function m(., x0), with x0 = 1/2 by default.

Subfolders include:

"BCCK": the method of Belloni et al. (2015).

"DW": the method of Dette and Wu (2026).

"New": the proposed method.

Run "marginal_mean/BCCK/main.R", "marginal_mean/DW/main.R", and "marginal_mean/New/Main.R" for the experiments associated with Tables S.16-S.17.  

5. Computational Efficiency (Supplementary Section S.4.6)

The folder "time_consuming" compares the running time of the proposed B-spline method with that of Dette and Wu (2026).

Run "time_consuming/main.R" for the comparison associated with Figure S.8.

The default settings are N = 100, T = 200, 500, 800, 1200, 1500, 2000, and five Monte Carlo replications for each sample size. The script stores the average timings in the R object "A" and displays a plot of the natural logarithm of running time.  

6. Time-Invariance Test (Supplementary Section S.10)

The folder "time_invariance" contains the small-sample size experiment for the corrected time-invariance test associated with Table S.18.

Run "time_invariance/main.R".

The test uses the centered contrast m_hat(tau, x) - m_hat_0(x), accounting for estimation of the unknown time-invariant mean. The default setting is (T, N) = (50, 400), with 500 Monte Carlo replications.

The output reports empirical type-I error rates for the original test and for the test with critical values multiplied by (1 + 20/T).  

7. Real Data Analyses

(1) Sydney Temperature

The folder "Sydney_temperature" contains the analysis in Section 7 of the main text, including Figures 1-2 and the corrected time-invariance test.

Run "Sydney_temperature/main.R". The input file "temperature_data.csv" is included and contains the data used for the 50 annual curves from 1959 to 2008.

The script produces the proposed 95% SCR with the finite-sample adjustment, computes the corrected time-invariance test, and calls "Wu.R" for the kernel-based comparison. It exports the estimated mean, upper surface, and lower surface as "z1.csv", "z2.csv", and "z3.csv". The corresponding kernel outputs are "wu1.csv", "wu2.csv", and "wu3.csv".

After running the R script, use "Sydney_temperature/plot.ipynb" to generate the surface plots from these CSV files. The notebook saves "sydney_temperature1.pdf" and "sydney_temperature_wu.pdf". 

For the finite-sample adjusted test, compare "test_stat" with "adjustfactor * Q_99" at the 1% level. 

(2) Australian Fertility Rates

The folder "Fertility" contains the analysis in Section S.3.1 of the Supplementary Material, including Figures S.1-S.4 and the additivity test.

Run "Fertility/main.R". The data are loaded from "Australiafertility" in the R package "rainbow" and cover 1921-2015 and maternal ages 15-49.

The script includes the bivariate mean estimator and 95% SCR, SCBs at fixed ages or years, and the additivity test. To obtain all marginal panels, rerun the fixed-age section with x0 = 25, 30, 35 and the fixed-year section with tau0 = 1930, 1960, 1990. The supplied selections are age 35 and year 1960. 

(3) Montreal Temperature

The folder "Montreal_temperature" contains the analysis in Section S.3.2 of the Supplementary Material, including Figures S.5-S.6 and the corrected time-invariance test.

Run "Montreal_temperature/main.R". The data are loaded from "MontrealTemp" in the R package "fda" and contain 34 annual curves from 1961 to 1994.

The script constructs the finite-sample adjusted 95% SCR, reports the corrected test statistic and critical values, and exports "z1.csv", "z2.csv", and "z3.csv". Run "Montreal_temperature/plot.ipynb" from the same folder to produce "Montreal-temp-95SCR.pdf".


Running environment

The analyses are written in R. Most entry scripts set the working directory through "rstudioapi::getActiveDocumentContext()". Open the desired entry script as the active document in RStudio before running it. Outside RStudio, replace that line with an explicit "setwd()" pointing to the folder containing the script.

The simulation entry scripts load local helper files, including "Gen_determinData.R", "Gen_randomdata.R", "test.R", "estimate.R", and the method-specific SCR/SCB functions. Keep the folder structure intact; an entry script is not a standalone file. Run each experiment in its own folder and preferably in a fresh R session, since different folders reuse function names and many drivers clear the workspace.

Ensure that the local "result" directory exists before running a simulation. If needed, run the following after setting the working directory:

    dir.create("result", showWarnings = FALSE)

Required packages

Packages used across the scripts include "Rcpp", "RcppEigen", "stringr", "tidyverse", "GeneralizedHyperbolic", "MASS", "matrixStats", "foreach", "doParallel", "zoo", "scatterplot3d", "plotly", "rgl", "rainbow", "fda", "ggplot2", and "rstudioapi". The alternative parallel drivers also use "pbapply". The R packages "splines" and "parallel" are also used.

The optional Python notebooks use "numpy", "pandas", and "matplotlib" and require a Jupyter-compatible environment. Run each notebook with its working directory set to the corresponding data folder. The notebooks contain "ax.w_xaxis", "ax.w_yaxis", and "ax.w_zaxis" calls; plotting compatibility depends on the installed Matplotlib version.

Notes

In the R scripts, "n" denotes the number of trajectories T, "N" denotes the within-curve grid size, and "MCMC" denotes the number of Monte Carlo datasets. 

If you encounter any issues, ensure that the mentioned R packages are installed. 

A random seed is specified in each script to reproduce the manuscript and supplemental results; 
slight variations may occur due to hardware differences, and the presented results have been reorganized from the computed outputs.

Simulation time depends on hardware; on a system with an 11th Gen Intel Core i7 CPU (2.80 GHz) and 16 GB RAM, a full run takes more than 24 hours. It's normal for the code to run without immediate feedback.

