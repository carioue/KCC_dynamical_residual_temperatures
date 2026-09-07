# KCC_dynamical_residual_temperatures

## Python scripts for training and fine-tuning
- Training_CMIP6_UV850_piCTL_all.py : pre-training on CMIP6 data (piControl simulations) with the wind at 850 hPa as predictor;
- Fine_Tuning_CMIP6_1850_1899.py : fine-tuning of the pretrained UNET on 1850-1899 of each model data;

## UNET trained on CMIP6 pre-industrial simualtions :
- Training_CMIP6_piCTL_all_UV850 : UNET pre-trained on CMIP6.

## R scripts for KCC
- KCC_Tdyn_Tres_WEU.R : code used to compute the observational constraints with KCC on dynamical and residual temperatures. See more details on this repository for instructions on how to use the KCC package https://gitlab.com/saidqasmi/KCC ;
- KCC_T_WEU.R : same but applied on T;
- KCC_files : Data used to run KCC scripts on T, $T_{Dyn}$ and $T_{Res}$;
- KCC_Tdyn_WEU_noNAT.nc, KCC_Tres_WEU.nc, KCC_T_WEU.nc : results of the constraints.
