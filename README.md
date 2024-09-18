# Spectral Deconfounding for High-Dimensional Sparse Additive Models

This repository contains the code for reproducing the plots in the paper
<i>Cyrill Scheidegger, Zijian Guo and Peter Bühlmann (2023). Spectral deconfounding for high-dimensional sparse additive models, arXiv:2312.02860</i>

The file `FitDeconfoundedHDAM.R`in the folder `FunctionsHDAM` contains the function `FitDeconfoundedHDAM` to fit a deconfounded high-dimensional additive model and some helper functions. The file `EstimationFactors.R` in the folder `FunctionsHDAM` contains the function `FitHDAM.withEstFactors` that implements an ad hoc method to achieve the same goal. The file `AnalyzeFittedHDAM.R` in the folder `FunctionsHDAM` contains functions to analyze/predict high-dimensional additive models based on the output of the function `FitDeconfoundedHDAM`.

The files in the folder `SimulationScripts` generate the plots from the paper. `ExampleIntroduction.R` generates Figure 1; `SingularValuePlot.R` generates Figure 2; `VarN.R` generates Figures 3, 4, 13 and 14; `VarP.R` generates Figures 5, 6, 15 and 16; `VarCS.R` generates Figures 7 and 8; `MotifEvaluation.R` generates Figures 9, 10, 11 and 12; `VarCP.R` generates Figures 17 and 18; `Nonlinear.R` generates Figures 19 and 20.
