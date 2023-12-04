# Spectral Deconfounding for High-Dimensional Sparse Additive Models

This repository contains the code for reproducing the plots in the paper "Cyrill Scheidegger, Zijian Guo and Peter Bühlmann (2023). Spectral deconfounding for high-dimensional sparse additive models"

The file `FitDeconfoundedHDAM.R`in the folder `FunctionsHDAM` contains the function `FitDeconfoundedHDAM` to fit a deconfounded high-dimensional additive model and some helper functions. The file `AnalyzeFittedHDAM.R` in the folder `FunctionsHDAM` contains functions to analyze/predict high-dimensional additive models based on the output of the function `FitDeconfoundedHDAM`.

The files in the folder `SimulationScripts` generate the plots from the paper. `ExampleIntroduction.R` generates Figure 1; `VarN.R`generates Figures 2 and 3; `VarP.R` generates Figures 4 and 5; `VaryCS` generates Figure 6; `VaryCP.R` generates Figure 7; `Nonlinear.R` generates Figure 8; `MotifEvaluation.R` generates Figures 9-12.
