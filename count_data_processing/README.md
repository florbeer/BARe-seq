This folder contains R scripts used for inference of the bursting parameters.
It is organized as follows:
- `R01_Optimization_Functions.R`: Script with custom function to infer the parameters for all the libraries. Functions were adjusted for the format of the count files for automated model fitting and parameter inference.
- `STAP / STARR ... NBfit.R` files: scripts that contain the parameter inference steps for each of the corresponding datasets, either promoter or enhancer libraries, as indicated by their name.
- `STAP / STARR ... mean_thres.R` files: scripts where mean expression threshold is determined for each library. The threshold is the point with minimal density when exploring the bimodal distribution of mean expression values.
- `SimulatedData_&_NBfit.R`: data simulation and parameter inference.
