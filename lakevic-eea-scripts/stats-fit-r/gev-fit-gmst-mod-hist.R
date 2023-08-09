library(extRemes)

setwd("C:\\Users\\rpietroi\\OneDrive - Vrije Universiteit Brussel\\repos_cloud\\lakevic-eea\\lakevic-eea-scripts")
getwd()

# Fitting models 

# 1. GMST only
# =============

GCMs = c('CanESM5','CNRM-CM6-1','GFDL-ESM4','IPSL-CM6A-LR','MIROC6','MRI-ESM2-0')

# test with just one model
#GCMs = GCMs[1]

# get fitted return period of observed intensity in 2020 from observations different possible

#a) 1897-2020 data (63.2 years) I think this is best one, most robust to excluding 2020 and downscaling temporal resolution
#rp_obs <- read.csv(paste("output/stat-fits-output/obs/lakevic-fit-obs-gmst-1897-2020.csv"), sep = ",", header = T)[1,"rp_2020"] 
#b) 1949-2020 data (48.5 years) I think this is less good 
rp_obs <- read.csv(paste("output/stat-fits-output/obs/lakevic-fit-obs-gmst-1949-2020.csv"), sep = ",", header = T)[1,"rp_2020"] 



for (i in GCMs){
  print(i)
  
  # get data
  dLdt <- read.csv(paste("output/event_def_dt180/hist-rcp70/blockmax_1850_2020_dt180_hist_", i, ".csv", sep=""), sep = ",", comment.char = "#", header = T, col.names = c("year", "dLdt", "NULL"))[,1:2]
  row.names(dLdt) <- dLdt$year
  gmst <- read.csv(paste("data/data-modified/gmst-models/hist-rcp370/wrt1850-1900/gmst_", i, "_hist-rcp370_1850_2100_wrt1850-1900_4yrlo.csv", sep="" ), sep = ",", header = T, col.names = c("year", "gmst"))[,1:2]
  row.names(gmst) <- gmst$year
  df <- merge(dLdt, gmst) # add iod later if necessary
  row.names(df) <- df$year
  
  # cut years
  startyear <- 1851
  endyear <- 2020
  slice <- as.character(c(startyear:endyear)) 
  df <- df[slice,]
  
  # fit GEV with only GMST as covariate (should give same output as Climate Explorer - indeed, it does)
  fit_gmst <- fevd(dLdt, df, location.fun = ~ gmst, type = "GEV", method = "MLE")
  
  # get matrix of covariate values for current & historic climates
  qcov_gmst_2020 <- make.qcov(fit_gmst, vals = list("mu1" = gmst["2020", "gmst"]))
  qcov_gmst_1900 <- make.qcov(fit_gmst, vals = list("mu1" = gmst["1900", "gmst"])) # or do -1.2
  
  # Save: model parameters
  params <- fit_gmst$results$par
  
  # nonstationary location & scale at given covariate value (note. location is sum)
  # Save loc in 2020 and 1900
  loc_2020 <- sum(findpars(fit_gmst, qcov = qcov_gmst_2020)$location)
  loc_1900 <- sum(findpars(fit_gmst, qcov = qcov_gmst_1900)$location)
  
  # calculate the return level in current & historic climate & hence change in intensity
  int_2020 <- return.level(fit_gmst, rp_obs, qcov = qcov_gmst_2020)[1]      # note: this is the fitted value, not the actual observation - should I use this or the actual obs?
  int_1900 <- return.level(fit_gmst, rp_obs, qcov = qcov_gmst_1900)[1]
  int_change <- (int_2020 - int_1900) 
  
  # calculate the return period of event in current & historic climate & hence probability ratio
  rp_2020 <- rp_obs                       # fixed by bias correction design
  rp_1900 <- 1/pextRemes(fit_gmst, int_2020, qcov = qcov_gmst_1900, lower.tail = F)
  prob_ratio <-  rp_1900 / rp_2020 
  
  # bootstrap probability ratio & change in intensity
  
  bootsize <-1000
  set.seed(1)
  boot_res <- sapply(1:bootsize, function(i) {                        # this sets the number of bootstrap samples
    
    boot_df <- df[sample(nrow(df), replace = T),]          # resample rows of dataframe (essentially resampling years)
    boot_mdl <- suppressWarnings(update(fit_gmst, data = boot_df))
    
    # calculate parameter CIs
    boot_params <- boot_mdl$results$par
    
    # calculate mu_2020 and mu_1900 CIs
    boot_loc_2020 <- sum(findpars(boot_mdl, qcov = qcov_gmst_2020)$location)
    boot_loc_1900 <- sum(findpars(boot_mdl, qcov = qcov_gmst_1900)$location)
    
    # calc change in intensity holding rp fixed to best estimate from observations
    boot_int_2020 <- return.level(boot_mdl, rp_2020, qcov = qcov_gmst_2020)[1]      # note: this is the fitted value, not the actual observation
    boot_int_1900 <- return.level(boot_mdl, rp_2020, qcov = qcov_gmst_1900)[1]
    boot_ic <- (boot_int_2020 - boot_int_1900)
    
    # calc return periods holding intensity fixed to best estimate from observations (incl. model bias correction)
    boot_rp_2020 <- unname(1/pextRemes(boot_mdl, int_2020, qcov = qcov_gmst_2020, lower.tail = F))
    boot_rp_1900 <- unname(1/pextRemes(boot_mdl, int_2020, qcov = qcov_gmst_1900, lower.tail = F))
    
    # calculate pr holding intensity fixed to best estimate from observations (incl. model bias correction)
    boot_pr <- unname(pextRemes(boot_mdl, int_2020, qcov = qcov_gmst_2020, lower.tail = F) / pextRemes(fit_gmst, int_2020, qcov = qcov_gmst_1900, lower.tail = F))
    
    
    c("params" = boot_params, "loc_2020" = boot_loc_2020, "loc_1900" = boot_loc_1900, 
      "rp_2020" = boot_rp_2020, "rp_1900" = boot_rp_1900 ,
      "PR" = boot_pr, "delta_I" = boot_ic, "int_2020" = boot_int_2020)
  })
  # find central 95% interval
  boot_ci <- apply(boot_res, 1, quantile, c(0.025, 0.975), na.rm = T)
  
  # make best estimate dataframe
  best_est <- c("params" = params, "loc_2020" = loc_2020, "loc_1900" = loc_1900, 
                "rp_2020" = unname(rp_2020), "rp_1900" = unname(rp_1900) ,
                "PR" = unname(prob_ratio), "delta_I" = int_change, "int_2020" = int_2020) 
  
  # add gmst vals and threshold from data in 2020
  gmst_2020 <- c(gmst["2020", "gmst"], NA, NA)
  gmst_1900 <- c(gmst["1900", "gmst"], NA, NA)
  
  # merge together
  res <- cbind(rbind(t(data.frame(best_est)), boot_ci), gmst_2020, gmst_1900)
  print(res)

  #save
  # write.csv(res, paste( "output/stat-fits-output/hist/rp",as.character(round(rp_obs,1)),"/lakevic-fit-hist-gmst-", as.character(i), "-rp", as.character(round(rp_obs,1)), "-", as.character(startyear), "-", as.character(endyear), ".csv", sep="" ), row.names=TRUE)
  
  
} 

