library(extRemes)

setwd("C:\\Users\\rpietroi\\OneDrive - Vrije Universiteit Brussel\\repos_cloud\\lakevic-eea\\lakevic-eea-scripts")
getwd()

# Fitting models - hist rp_2020 vs hist-nat rp_2020 

# 1. GMST only
# =============

GCMs = c('CanESM5','CNRM-CM6-1','GFDL-ESM4','IPSL-CM6A-LR','MIROC6','MRI-ESM2-0')

# test with just one model
#GCMs = GCMs[1]
#i <- GCMs


# get fitted return period of observed intensity in 2020 from observations different possible

#a) 1897-2020 data (63.2 years)
#rp_obs <- read.csv(paste("output/stat-fits-output/obs/lakevic-fit-obs-gmst-1897-2020.csv"), sep = ",", header = T)[1,"rp_2020"] 
#b) 1949-2020 data (48.5 years) 
rp_obs <- read.csv(paste("output/stat-fits-output/obs/lakevic-fit-obs-gmst-1949-2020.csv"), sep = ",", header = T)[1,"rp_2020"] 



for (i in GCMs){
  print(i)
  
  # get data
  dLdt_h <- read.csv(paste("output/event_def_dt180/hist-rcp70/blockmax_1850_2020_dt180_hist_", i, ".csv", sep=""), sep = ",", comment.char = "#", header = T, 
                     col.names = c("year", "dLdt_h", "NULL"))[,1:2]
  row.names(dLdt_h) <- dLdt_h$year
  dLdt_hn <- read.csv(paste("output/event_def_dt180/hist-nat/blockmax_1850_2020_dt180_hist-nat_", i, ".csv", sep=""), sep = ",", comment.char = "#", header = T, 
                      col.names = c("year", "dLdt_hn", "NULL"))[,1:2]
  row.names(dLdt_hn) <- dLdt_hn$year
  gmst_h <- read.csv(paste("data/data-modified/gmst-models/hist-rcp370/wrt1850-1900/gmst_", i, "_hist-rcp370_1850_2100_wrt1850-1900_4yrlo.csv", sep="" ), sep = ",", header = T, 
                     col.names = c("year", "gmst_h"))[,1:2]
  row.names(gmst_h) <- gmst_h$year
  gmst_hn <- read.csv(paste("data/data-modified/gmst-models/hist-nat/wrt1850-1900/gmst_", i, "_hist-nat_1850_2020_wrt1850-1900_4yrlo.csv", sep="" ), sep = ",", header = T, 
                      col.names = c("year", "gmst_hn"))[,1:2]
  row.names(gmst_hn) <- gmst_hn$year
  df <- merge( merge( merge(dLdt_h, dLdt_hn ), gmst_h), gmst_hn) # add iod later if necessary
  row.names(df) <- df$year
  
  # cut years
  startyear <- 1851
  endyear <- 2020
  slice <- as.character(c(startyear:endyear)) 
  df <- df[slice,]
  
  
  # a) hist: fit GEV with only GMST as covariate to hist
  fit_gmst_h <- fevd(dLdt_h, df, location.fun = ~ gmst_h, type = "GEV", method = "MLE")
  
  # get matrix of covariate values for current climates (hist simulation)
  qcov_gmst_2020_h <- make.qcov(fit_gmst_h, vals = list("mu1" = gmst_h["2020", "gmst_h"]))
  
  # Check: model parameters are the same as shift fit to hist alone
  params_h <- fit_gmst_h$results$par
  
  # calc return level in current climate in hist (model-specific intensity threshold for rp_obs)
  int_mod <- return.level(fit_gmst_h, rp_obs, qcov = qcov_gmst_2020_h)[1]
  int_2020_h <- int_mod
  
  # return period in hist fit is by definition rp_obs
  rp_2020_h <- unname(rp_obs) # 1/pextRemes(fit_gmst_h, int_mod, qcov = qcov_gmst_2020_h, lower.tail = F)
  
  
  # b) hist-nat: calculate return period of model-specific intensity threshold in hist-nat fit 
  fit_gmst_hn <- fevd(dLdt_hn, df, location.fun = ~ gmst_hn, type = "GEV", method = "MLE")
  
  # Check: model parameters 
  params_hn <- fit_gmst_hn$results$par
  
  # get matrix of covariate values for current climates
  qcov_gmst_2020_hn <- make.qcov(fit_gmst_hn, vals = list("mu1" = gmst_hn["2020", "gmst_hn"]))

  # calculate the return period of event (int_mod) in 2020 climate (hist-nat) and probability ratio
  rp_2020_hn <- unname(1/pextRemes(fit_gmst_hn, int_mod, qcov = qcov_gmst_2020_hn, lower.tail = F)) 
  prob_ratio <- rp_2020_hn / rp_2020_h
  
  # calculate return level of event (rp_obs) in current climate (hist-nat) and change in intensity
  int_2020_hn <- return.level(fit_gmst_hn, rp_obs, qcov = qcov_gmst_2020_hn)[1]
  int_change <- int_2020_h - int_2020_hn
  
  
  # c) bootstrap PR and change in intensity
  bootsize <- 1000
  set.seed(1)
  boot_res <- sapply(1:bootsize, function(i) {
    
    boot_df <- df[sample(nrow(df), replace = T),]          # resample rows of dataframe (essentially resampling years)
    boot_mdl_h <- suppressWarnings(update(fit_gmst_h, data = boot_df))
    boot_mdl_hn <- suppressWarnings(update(fit_gmst_hn, data = boot_df))
    
    # calc return periods holding intensity fixed and pr 
    boot_rp_2020_h <- unname(1/pextRemes(boot_mdl_h, int_mod, qcov = qcov_gmst_2020_h, lower.tail = F))
    boot_rp_2020_hn <-  unname(1/pextRemes(boot_mdl_hn, int_mod, qcov = qcov_gmst_2020_hn, lower.tail = F))
    boot_pr <- boot_rp_2020_hn / boot_rp_2020_h
    
    # calc return level holding rp fixed to rp_obs (this is correct, right?)
    boot_int_2020_h <- return.level(boot_mdl_h, rp_obs, qcov = qcov_gmst_2020_h)[1] 
    boot_int_2020_hn <- return.level(boot_mdl_hn, rp_obs, qcov = qcov_gmst_2020_hn)[1] 
    boot_ic <- (boot_int_2020_h - boot_int_2020_hn)
    
    c("rp_2020_h" = boot_rp_2020_h, "rp_2020_hn" = boot_rp_2020_hn,
      "int_2020_h" = boot_int_2020_h, "int_2020_hn" = boot_int_2020_hn,
      "PR" = boot_pr, "delta_I" = boot_ic)
    
  })
  
  # find central 95% interval
  boot_ci <- apply(boot_res, 1, quantile, c(0.025, 0.975), na.rm = T)
  
  # make best estimate dataframe
  best_est <- c("rp_2020_h" = rp_2020_h, "rp_2020_hn" = rp_2020_hn,
                "int_2020_h" = int_2020_h, "int_2020_hn" = int_2020_hn,
                "PR" = prob_ratio, "delta_I" = int_change)
  # merge together
  res <- rbind(t(data.frame(best_est)), boot_ci)
  print(res)
  
  #save
  write.csv(res, paste( "output/stat-fits-output/hist-vs-hist-nat/rp",as.character(round(rp_obs,1)),"/lakevic-fit-hist-vs-histnat-gmst-", as.character(i), "-rp", as.character(round(rp_obs,1)), "-", as.character(startyear), "-", as.character(endyear), ".csv", sep="" ), row.names=TRUE)
  
  
}

