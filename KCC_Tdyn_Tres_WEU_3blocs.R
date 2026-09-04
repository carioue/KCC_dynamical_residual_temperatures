# Loading KCC package
library(KCC,lib.loc='/KCC/packages/4.3')
# and other useful package(s)
library(abind,lib.loc='/KCC/packages/4.3')
library(ncdf4)


set.seed(1)

do_antnat_dec_glo = F
do_antnat_dec_loc = F 
do_fit_MAR_obs_glo = F
do_fit_MAR_obs_loc = F

# Sample size to derive normal distributions
Nres = 1000
sample_str = c("be", paste0("nres",1:Nres))

year_g = 1850:2100
ny_g = length(year_g)

load("/KCC/packages/FF_CMIP6.rda")

# Load EBM parameters fitted on available CMIP6 models
load("/KCC/packages/ebm_params.rda")

#ebm_params

# Compute the EBM response e, to be used in equation 7 in Qasmi and Ribes (2021)
e = ebm_response(FF,ebm_params,year_g,Nres)

# GLO MODELS
if (do_antnat_dec_glo == T) {

	# Load historical+ssp585 GMST ensemble means from CMIP6 models
	X_glo_scen_tmp = loadRData("Tas_glo_CMIP6_histssp585_ann.Rdata")

	# Extract the scenario and the chosen simulated years
	X_glo_scen_full = X_glo_scen_tmp[as.character(year_g),,]

	# Consider the same models in historical+ssp585 and piControl simulations
	models_scen_full = dimnames(X_glo_scen_full)$model
	# Check if any NA in times series
	models = models_scen_full[apply(is.na(X_glo_scen_full[,,"value"]), 2, sum) == 0] 
	# The following CMIP6 models will be used to derive the prior
	models

	Nmod = length(models)
	X_glo_scen = X_glo_scen_full[,models,]

	X_fit_glo = array(NA, dim = c(ny_g, 2, Nmod),
		  dimnames = list(year = as.character(year_g),
		                forcing = c("all","nat"),
		                model = models))
	# degrees of freedom for the spline function
	df = 6                            
	message("Decomposition for each model (take some time!)...")
	X_fit_glo[,c("all","nat"),] = x_fit(X_glo_scen[,,"value"], e, df, ant=F)
	save("X_fit_glo",file="X_fit_glo_ann.Rdata")
	
} else {

	load("X_fit_glo_ann.Rdata")
}
#GLO OBS
# Load observed GMST from Hadcrut over 1850-2024
year_obs_glo = 1850:2024
ref_glo = 1961:1990 
ny_obs_glo = length(year_obs_glo)
Xo_glo_full = loadRData("Xo_raw.Rdata")
Xo_glo = Xo_glo_full[as.character(year_obs_glo),]

if (do_fit_MAR_obs_glo == T) {
	 # Estimate the response to all external forcings by the multimodel ensemble mean
	raw_mmm = apply(X_fit_glo[,"all",], 1, mean, na.rm=T)
	# Calculate obs residuals, our estimate of internal variability
	raw_mmm_c = raw_mmm - mean(raw_mmm[year_g %in% ref_glo]) # raw_mmm_c must be anomalies wrt 1961-1990
	Xo_glo_med = apply(Xo_glo, 1, median)
	Xo_glo_c = Xo_glo_med - mean(Xo_glo_med[year_obs_glo %in% ref_glo])
	res_glo = Xo_glo_c - raw_mmm_c[as.character(year_obs_glo)]
	# Fit the parameters of the MAR models on residuals
	message("Fitting MAR parameters (may take some time!)...")
	theta_obs_glo = estim_mar2_link(res_glo)
	# Compute the associated covariance matrix
	Sigma2_obs_iv_glo = Sigma_mar2(theta_obs_glo,res_glo)
		
	# Check that MAR modelling is ok
	acf_res = acf(res_glo, plot=F)
	par(cex=1.1,font=2, font.lab=2,font.axis=2,lwd=2,mgp=c(2.8,.7,0),mar=c(4,4,1,4),las=2,tcl=-.4,cex.lab=1.2)
	plot(acf_res,ylab="ACF of GMST residuals wrt CMIP6 multimodel mean",ci=0,lwd=2)
	colors_acf = c("red", "black")
	lines(0:20, Sigma2_obs_iv_glo[1,1:21]/Sigma2_obs_iv_glo[1,1], col=colors_acf[1], lwd=2)
	legend("topright", legend = c("MAR fit", "OBS"), col = colors_acf, lwd=2, lty=1)

	# Add the covariance matrix for measurement errors
	Sigma2_obs_glo = Sigma2_obs_iv_glo + var(t(Xo_glo))
	
	save(res_glo, theta_obs_glo, Sigma2_obs_glo, file="MAR_obs_glo_ann.Rdata")
	
} else {

	load("MAR_obs_glo_ann.Rdata")

}

year_l = 1900:2100
ny_l = length(year_l)
era5_dyn = nc_open('Tdyn_era5_day_1940_2024_yearmean_fldmean_WEU_land.nc')
era5_res = nc_open('Tres_era5_day_1940_2024_yearmean_fldmean_WEU_land.nc')
Xo_loc_full_dyn = ncvar_get(era5_dyn,"Tdyn")
Xo_loc_full_res = ncvar_get(era5_res,"Tres")
year_obs_loc = 1940:2024
ny_obs_loc = length(year_obs_loc)
ref_loc=1940:2024

e_loc = ebm_response(FF,ebm_params,year_l,Nres)

#LOC MODELS
if (do_antnat_dec_loc == T) {
	# Load historical+ssp585 GMST ensemble means from CMIP6 models
	X_loc_scen_dyn_tmp = loadRData("Tdyn_WEU_CMIP6_histssp585_ann.Rdata")
	X_loc_scen_res_tmp = loadRData("Tres_WEU_CMIP6_histssp585_ann.Rdata")
	# Extract the scenario and the chosen simulated years
	X_loc_scen_full_dyn = X_loc_scen_dyn_tmp[as.character(year_l),,]
	X_loc_scen_full_res = X_loc_scen_res_tmp[as.character(year_l),,]

	# Consider the same models in historical+ssp585 and piControl simulations
	models_scen_full = dimnames(X_loc_scen_full_dyn)$model
	# Check if any NA in times series
	models = models_scen_full[apply(is.na(X_loc_scen_full_dyn[,,"value"]), 2, sum) == 0] 
	# The following CMIP6 models will be used to derive the prior
	models

	Nmod = length(models)
	X_loc_scen_dyn = X_loc_scen_full_dyn[,models,]
	X_loc_scen_res = X_loc_scen_full_res[,models,]

	X_fit_loc_dyn = array(NA, dim = c(ny_l, 2, Nmod),
		  dimnames = list(year = as.character(year_l),
				forcing = c("all","nat"),
				model = models))
	X_fit_loc_res=X_fit_loc_dyn
	# degrees of freedom for the spline function
	df = 6                            
	message("Decomposition for each model (take some time!)...")
	X_fit_loc_dyn[,c("all","nat"),] = x_fit(X_loc_scen_dyn[,,"value"], e_loc, df, ant=F)
	X_fit_loc_res[,c("all","nat"),] = x_fit(X_loc_scen_res[,,"value"], e_loc, df, ant=F)
	save("X_fit_loc_dyn",file=("X_fit_loc_dyn.Rdata"))
	save("X_fit_loc_res",file=("X_fit_loc_res.Rdata"))
	
} else {

	load(paste0("X_fit_loc_dyn.Rdata"))
	load(paste0("X_fit_loc_res.Rdata"))

}

#LOC OBS
Xo_loc_dyn = Xo_loc_full_dyn
Xo_loc_res = Xo_loc_full_res
names(Xo_loc_dyn) = year_obs_loc
names(Xo_loc_res) = year_obs_loc

if (do_fit_MAR_obs_loc == T) {
	 # Estimate the response to all external forcings by the multimodel ensemble mean
	raw_mmm = apply(X_fit_loc_dyn[,"all",], 1, mean, na.rm=T)
	# Calculate obs residuals, our estimate of internal variability
	raw_mmm_c = raw_mmm - mean(raw_mmm[year_l %in% ref_loc])
	Xo_loc_med = Xo_loc_dyn
	Xo_loc_c = Xo_loc_med - mean(Xo_loc_med[year_obs_loc %in% ref_loc])
	res_loc = Xo_loc_c - raw_mmm_c[as.character(year_obs_loc)]
	# Fit the parameters of the MAR models on residuals
	message("Fitting MAR parameters (may take some time!)...")
	theta_obs_loc = estim_mar2_link(res_loc)
	# Compute the associated covariance matrix
	Sigma2_obs_iv_loc_dyn = Sigma_mar2(theta_obs_loc,res_loc)
		
	# Check that MAR modelling is ok
	acf_res = acf(res_loc, plot=F)
	par(cex=1.1,font=2, font.lab=2,font.axis=2,lwd=2,mgp=c(2.8,.7,0),mar=c(4,4,1,4),las=2,tcl=-.4,cex.lab=1.2)
	plot(acf_res,ylab="ACF of GMST residuals wrt CMIP6 multimodel mean",ci=0,lwd=2)
	colors_acf = c("red", "black")
	lines(0:20, Sigma2_obs_iv_loc_dyn[1,1:21]/Sigma2_obs_iv_loc_dyn[1,1], col=colors_acf[1], lwd=2)
	legend("topright", legend = c("MAR fit", "OBS"), col = colors_acf, lwd=2, lty=1)

	# Add the covariance matrix for measurement errors
	Sigma2_obs_loc_dyn = Sigma2_obs_iv_loc_dyn #+ var(t(Xo_loc))
	save(res_loc, theta_obs_loc, Sigma2_obs_iv_loc_dyn, file=paste0("MAR_obs_tas_loc_dyn.Rdata"))
	
	# Estimate the response to all external forcings by the multimodel ensemble mean
	raw_mmm = apply(X_fit_loc_res[,"all",], 1, mean, na.rm=T)
	# Calculate obs residuals, our estimate of internal variability
	raw_mmm_c = raw_mmm - mean(raw_mmm[year_l %in% ref_loc])
	Xo_loc_med = Xo_loc_res
	Xo_loc_c = Xo_loc_med - mean(Xo_loc_med[year_obs_loc %in% ref_loc])
	res_loc = Xo_loc_c - raw_mmm_c[as.character(year_obs_loc)]
	# Fit the parameters of the MAR models on residuals
	message("Fitting MAR parameters (may take some time!)...")
	theta_obs_loc = estim_mar2_link(res_loc)
	# Compute the associated covariance matrix
	Sigma2_obs_iv_loc_res = Sigma_mar2(theta_obs_loc,res_loc)
		
	# Check that MAR modelling is ok
	acf_res = acf(res_loc, plot=F)
	par(cex=1.1,font=2, font.lab=2,font.axis=2,lwd=2,mgp=c(2.8,.7,0),mar=c(4,4,1,4),las=2,tcl=-.4,cex.lab=1.2)
	plot(acf_res,ylab="ACF of GMST residuals wrt CMIP6 multimodel mean",ci=0,lwd=2)
	colors_acf = c("red", "black")
	lines(0:20, Sigma2_obs_iv_loc_res[1,1:21]/Sigma2_obs_iv_loc_res[1,1], col=colors_acf[1], lwd=2)
	legend("topright", legend = c("MAR fit", "OBS"), col = colors_acf, lwd=2, lty=1)

	# Add the covariance matrix for measurement errors
	Sigma2_obs_loc_res = Sigma2_obs_iv_loc_res #+ var(t(Xo_loc))
	
	save(res_loc, theta_obs_loc, Sigma2_obs_iv_loc_res, file=paste0("MAR_obs_tas_loc_res.Rdata"))
	
} else {

	load(paste0("MAR_obs_tas_loc_dyn.Rdata"))
	load(paste0("MAR_obs_tas_loc_res.Rdata"))

}

X_fit = abind(X_fit_loc_dyn,X_fit_loc_res, X_fit_glo, along=1, use.dnns=T)
# Rearrange dimensions along year dimension
year_loc_dyn = paste0(as.character(year_l),"locdyn")
year_loc_res = paste0(as.character(year_l),"locres")
year_glo = paste0(as.character(year_g),"glo")
dimnames(X_fit)$year = c(year_loc_dyn,year_loc_res, year_glo)

# Bind arrays containing the local and global observations
Xo = abind(replicate(200,Xo_loc_dyn),replicate(200,Xo_loc_res), Xo_glo[,1:200], along=1, use.dnns=T)
year_obs_loc_dyn_str = paste0(as.character(year_obs_loc),"locdyn")
year_obs_loc_res_str = paste0(as.character(year_obs_loc),"locres")
year_obs_glo_str = paste0(as.character(year_obs_glo),"glo")
dimnames(Xo)$year = c(year_obs_loc_dyn_str,year_obs_loc_res_str, year_obs_glo_str)


# Compute the associated covariance matrix
#Sigma2_obs_iv_loc_glo = Sigma_mar_dep(theta_obs_loc_glo,res_loc_glo)
Sigma2_obs_iv_loc_glo = matrix(0,nrow=dim(Xo)[1],ncol=dim(Xo)[1])
Sigma2_obs_iv_loc_glo[1:ny_obs_loc,1:ny_obs_loc] = Sigma2_obs_iv_loc_dyn
Sigma2_obs_iv_loc_glo[ny_obs_loc+(1:ny_obs_loc),ny_obs_loc+(1:ny_obs_loc)] = Sigma2_obs_iv_loc_res
Sigma2_obs_iv_loc_glo[2*ny_obs_loc+(1:ny_obs_glo),2*ny_obs_loc+(1:ny_obs_glo)] = Sigma2_obs_iv_glo

Sigma2_obs_loc_glo = Sigma2_obs_iv_loc_glo + var(t(Xo))

# Constrain by local+global observations
X_krig_loc_glo_list = prior2posterior(X_fit, Xo, Sigma2_obs_loc_glo, Nres, centering_CX=T, ref_CX=c(paste0(1940:2024,"locdyn"),paste0(1940:2024,"locres"),paste0(1961:1990,"glo"))) 

# Convert output lists to arrays
# X_unconstrained
X_uncons_loc_glo = mvgauss_to_Xarray(X_krig_loc_glo_list$uncons$mean,X_krig_loc_glo_list$uncons$var,Nres)
# X_constrained
X_cons_loc_glo = mvgauss_to_Xarray(X_krig_loc_glo_list$cons$mean,X_krig_loc_glo_list$cons$var,Nres)

# Bind X_uncons and X_cons together using the abind library
l_tmp = c(dimnames(X_uncons_loc_glo)[1:3], list(constrain=c("uncons","cons")))
X_krig_loc_glo = abind(X_uncons_loc_glo, X_cons_loc_glo, along=4)
dimnames(X_krig_loc_glo) = l_tmp	
X_krig_tdyn = X_krig_loc_glo[1:201,,"all",]
X_krig_tres = X_krig_loc_glo[202:402,,"all",]
#stop('stop')

save(X_krig_tdyn, file=paste0("KCC_Tdyn_WEU.Rdata"))
save(X_krig_tres, file=paste0("KCC_Tres_WEU.Rdata"))


