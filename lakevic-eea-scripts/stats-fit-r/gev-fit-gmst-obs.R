library(extRemes)

setwd("C:\\Users\\rpietroi\\OneDrive - Vrije Universiteit Brussel\\repos_cloud\\lakevic-eea\\lakevic-eea-scripts")
getwd()

# Fitting models 

# 1. GMST only
# =============

# load data. Different obs lake level data sources : 
# blockmax_1897_2021_dt180_obs.csv
# blockmax_1948_2022_dt180_obs.csv (this and the previous are the same for overlapping years)
# blockmax_1897_2021_dt180_obs_downscaled.csv

#dLdt <- read.csv("output/event_def_dt180/obs/blockmax_1897_2021_dt180_obs.csv", sep = ",", comment.char = "#", header = T, col.names = c("year", "dLdt", "NULL"))[,1:2]
#dLdt <- read.csv("output/event_def_dt180/obs/blockmax_1897_2021_dt180_obs_downscaled.csv", sep = ",", comment.char = "#", header = T, col.names = c("year", "dLdt", "NULL"))[,1:2]
dLdt <- read.csv("output/event_def_dt180/obs/blockmax_1897_2021_dt180_obs_rm_overlap.csv", sep = ",", comment.char = "#", header = T, col.names = c("year", "dLdt", "NULL"))[,1:2]

row.names(dLdt) <- dLdt$year

gmst <- read.csv("data/input-data/gmst-obs/igiss_al_gl_a_4yrlo.csv", sep = ",", comment.char = "#", header = T, col.names = c("year", "gmst"))[,1:2]
row.names(gmst) <- gmst$year

iod <- read.csv("data/data-modified/iod-obs/dmi_ond_1870_2021.csv", sep = ",", comment.char = "#", header = T, col.names = c("year", "iod"))[,1:2]
row.names(iod) <- iod$year

df <- merge(dLdt, merge(gmst,iod))
row.names(df) <- df$year

# cut to desired years. Different possible: 1949-2019, 1949-2020, 1897-2020, 1897-2019
startyear <- 1897 # 1949 # 1897
endyear <- 2020 # 20 # 2019
slice <- as.character(c(startyear:endyear)) 
df <- df[slice,]

# fit GEV with only GMST as covariate (should give same output as Climate Explorer - indeed, it does)
fit_gmst <- fevd(dLdt, df, location.fun = ~ gmst, type = "GEV", method = "MLE")

# only for the dataset with missing data 
fit_gmst <- fevd(dLdt, df, location.fun = ~ gmst, type = "GEV", method = "MLE", 
                 na.action = na.exclude)

# get matrix of covariate values for current & historic climates
qcov_gmst_2020 <- make.qcov(fit_gmst, vals = list("mu1" = gmst["2020", "gmst"]))
qcov_gmst_1900 <- make.qcov(fit_gmst, vals = list("mu1" = gmst["1900", "gmst"])) # or do -1.2

# Save: model parameters
params <- fit_gmst$results$par

# nonstationary location & scale at given covariate value (note. location is sum)
findpars(fit_gmst, qcov = qcov_gmst_2020)
findpars(fit_gmst, qcov = qcov_gmst_1900)

# Save loc in 2020 and 1900
loc_2020 <- sum(findpars(fit_gmst, qcov = qcov_gmst_2020)$location)
loc_1900 <- sum(findpars(fit_gmst, qcov = qcov_gmst_1900)$location)

# calculate the return period of 2022 precip in current & historic climate & hence probability ratio
rp_2020 <- 1/pextRemes(fit_gmst, dLdt["2020", "dLdt"], qcov = qcov_gmst_2020, lower.tail = F)      # remember to tell it whether to look at the lower or upper tail!
rp_1900 <- 1/pextRemes(fit_gmst, dLdt["2020", "dLdt"], qcov = qcov_gmst_1900, lower.tail = F)
prob_ratio <-  rp_1900 / rp_2020  

# calculate the return level in current & historic climate & hence change in intensity
int_2020 <- return.level(fit_gmst, rp_2020, qcov = qcov_gmst_2020)[1]      # note: this is the fitted value, not the actual observation - should I use this or the actual obs?
int_1900 <- return.level(fit_gmst, rp_2020, qcov = qcov_gmst_1900)[1]
int_change <- (int_2020 - int_1900) 


# bootstrap probability ratio & change in intensity

bootsize <-10000 #1000 to test !

set.seed(1)
boot_res <- sapply(1:bootsize, function(i) {                        # this sets the number of bootstrap samples
  
  boot_df <- df[sample(nrow(df), replace = T),]          # resample rows of dataframe (essentially resampling years)
  boot_mdl <- suppressWarnings(update(fit_gmst, data = boot_df))
  
  # calculate parameter CIs
  boot_params <- boot_mdl$results$par
  
  # calculate mu_2020 and mu_1900 CIs
  boot_loc_2020 <- sum(findpars(boot_mdl, qcov = qcov_gmst_2020)$location)
  boot_loc_1900 <- sum(findpars(boot_mdl, qcov = qcov_gmst_1900)$location)
  
  # calc return periods
  boot_rp_2020 <- unname(1/pextRemes(boot_mdl, dLdt["2020", "dLdt"], qcov = qcov_gmst_2020, lower.tail = F))
  boot_rp_1900 <- unname(1/pextRemes(boot_mdl, dLdt["2020", "dLdt"], qcov = qcov_gmst_1900, lower.tail = F))
  
  # calculate pr & change in intensity
  boot_pr <- unname(pextRemes(boot_mdl, dLdt["2020", "dLdt"], qcov = qcov_gmst_2020, lower.tail = F) / pextRemes(fit_gmst, dLdt["2020", "dLdt"], qcov = qcov_gmst_1900, lower.tail = F))
  
  boot_int_2020 <- return.level(boot_mdl, rp_2020, qcov = qcov_gmst_2020)[1]      # note: this is the fitted value, not the actual observation
  boot_int_1900 <- return.level(boot_mdl, rp_2020, qcov = qcov_gmst_1900)[1]
  boot_ic <- (boot_int_2020 - boot_int_1900) 
  
  c("params" = boot_params, "loc_2020" = boot_loc_2020, "loc_1900" = boot_loc_1900, 
    "rp_2020" = boot_rp_2020, "rp_1900" = boot_rp_1900 ,
    "PR" = boot_pr, "delta_I" = boot_ic)
})
# find central 95% interval
boot_ci <- apply(boot_res, 1, quantile, c(0.025, 0.975), na.rm = T)
boot_ci

# make best estimate dataframe
best_est <- c("params" = params, "loc_2020" = loc_2020, "loc_1900" = loc_1900, 
                "rp_2020" = unname(rp_2020), "rp_1900" = unname(rp_1900) ,
                "PR" = unname(prob_ratio), "delta_I" = int_change) 

# add gmst vals and threshold from data in 2020
gmst_2020 <- c(gmst["2020", "gmst"], NA, NA)
gmst_1900 <- c(gmst["1900", "gmst"], NA, NA)
int_obs <- c(dLdt["2020", "dLdt"], NA, NA) # for models add modelled intensity instead of this, and calc also the CI

# merge together
res <- cbind(rbind(t(data.frame(best_est)), boot_ci), int_obs, gmst_2020, gmst_1900)

#save

#write.csv(res, paste( "output/stat-fits-output/lakevic-fit-obs-gmst-", as.character(startyear), "-", as.character(endyear), ".csv", sep="" ), row.names=TRUE)

#write.csv(res, paste( "output/stat-fits-output/lakevic-fit-obs-downsc-gmst-", as.character(startyear), "-", as.character(endyear), ".csv", sep="" ), row.names=TRUE)

#write.csv(res, paste( "output/stat-fits-output/lakevic-fit-obs-gmst-boot10k-", as.character(startyear), "-", as.character(endyear), ".csv", sep="" ), row.names=TRUE)

#write.csv(res, paste( "output/stat-fits-output/lakevic-fit-obs-rm-overlap-gmst-", as.character(startyear), "-", as.character(endyear), ".csv", sep="" ), row.names=TRUE)

