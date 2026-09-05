# KCC_dynamical_residual_temperatures

## Python scripts for training and fine-tuning
- Training_CMIP6_UV850_piCTL.py : pre-training on CMIP6 data (piControl simulations) with the wind at 850 hPa as predictor;
- FineTuning_ERA519401969_UV850.py : fine-tuning of the pretrained UNET (with UV850) on ERA5 1940-1969

All the UNETs fine-tuned on models can be found on the Zenodo link, as well as the dynamical and residual contributions estimated for all the simulations.
  
## R scripts for KCC
- KCC_Tdyn_Tres_WEU.R : code used to compute the observational constraints with KCC on dynamical and residual temperatures. See more details on this repository for instructions on how to use the KCC package https://gitlab.com/saidqasmi/KCC .
- KCC_files : Data used to run KCC script on dynamical and residual temperatures and the results of the constraints (for T, $T_{Dyn}$ and $T_{Res}$).
