This folder contains R scripts used for inference of the bursting parameters.
It is organized as follows:
- `R01_Optimization_Functions.R`: Script with custom function to infer the parameters for all the libraries. Functions were adjusted for the format of the count files for automated model fitting and parameter inference.
- `STAP / STARR ... NBfit.R` files: scripts that contain the parameter inference steps for each of the corresponding datasets, either promoter or enhancer libraries, as indicate dby their name
- `SimulatedData_&_NBfit.R`: data simulation and parameter inference
