library(ggplot2)
library(mgcv)
library(plyr)
library(dplyr)

### This model only requires annual survival estimates - accounting for transience in survival and trap dependence and with time-dependent recapture (see Supplementary information S6).
# SST effects are fitted on the logit transformed raw data and predicted under current and future SSTs 

rm(list = ls())

# 
# 1. Read survival estimates
# 
pdata <- read.csv("Data_Surv.csv")

# Fix encoding with norwegian characters
pdata[] <- lapply(pdata, function(x) {
  if (is.character(x)) iconv(x, from = "latin1", to = "UTF-8") else x
})

# Compute logit-scale response + weights
pdata$phi_logit  <- qlogis(pdata$mean)
pdata$phi_var    <- (pdata$se / (pdata$mean * (1 - pdata$mean)))^2
pdata$phi_weight <- 1 / pdata$phi_var

#
# 2. Read SST data (autumn and winter)
# 
covA <- read.csv("Data_SSTaut.csv", sep = ";")
covW <- read.csv("Data_SSTwin.csv", sep = ";")

# Fix encoding
covA[] <- lapply(covA, \(x) if(is.character(x)) iconv(x,"latin1","UTF-8") else x)
covW[] <- lapply(covW, \(x) if(is.character(x)) iconv(x,"latin1","UTF-8") else x)

covA$sp_col2 <- paste(covA$species, covA$site, sep = "-")
covW$sp_col2 <- paste(covW$species, covW$site, sep = "-")
pdata$sp_col2 <- paste(pdata$species, pdata$colony, sep = "-")

# Merge SST values
pdata$sst_aut_raw <- covA$MeanAut_SST[ match(pdata$sp_col2, covA$sp_col2) ]
pdata$sst_win_raw <- covW$MeanWin_SST[ match(pdata$sp_col2, covW$sp_col2) ]

# 
# 3. Scale SST within population
#
scale_by_pop <- function(x, grp) {
  x_scaled <- x
  attr_center <- numeric(length(unique(grp)))
  attr_scale  <- numeric(length(unique(grp)))
  
  u <- unique(grp)
  
  for (i in seq_along(u)) {
    idx <- grp == u[i]
    x_scaled[idx] <- scale(x[idx])
    attr_center[i] <- attr(scale(x[idx]), "scaled:center")
    attr_scale[i]  <- attr(scale(x[idx]), "scaled:scale")
  }
  
  list(x_scaled = x_scaled,
       centre = setNames(attr_center, u),
       scale  = setNames(attr_scale,  u))
}

# Autumn scaling
scA <- scale_by_pop(pdata$sst_aut_raw, pdata$sp_col2)
pdata$sst_aut <- scA$x_scaled

# Winter scaling
scW <- scale_by_pop(pdata$sst_win_raw, pdata$sp_col2)
pdata$sst_win <- scW$x_scaled


# -------------------------------------------------------------
# 4. Fit models separately for Autumn & Winter
# -------------------------------------------------------------

predict_survival <- function(model, beta_vcv, new_sst_scaled) {
  # linear predictor
  x <- coef(model)[1] + coef(model)[2] * new_sst_scaled
  
  # survival
  surv <- plogis(x)
  
  # delta-method SE
  mat <- c(1, new_sst_scaled)
  var <- t(mat) %*% beta_vcv %*% mat
  se <- sqrt(var)
  
  list(mu = surv, sd = se)
}

# -------------------------------------------------------------
# 5. For each population: compute CURRENT vs FUTURE predictions
# -------------------------------------------------------------
results_out <- data.frame()

pops <- unique(pdata$sp_col2)

for (pop in pops) {
  
  dsub <- pdata[pdata$sp_col2 == pop, ]
  sp  <- dsub$species[1]
  col <- dsub$colony[1]
  
  # Choose season based on available SST columns
  # Here we show AUTUMN example ??? duplicate for winter
  # -------------------------------------------------------------------
  # Fit model (AUTUMN)
  modA <- lm(phi_logit ~ sst_aut, data = dsub, weights = phi_weight)
  
  # Extract covariance matrix
  vcvA <- vcov(modA)
  
  # Scaling used
  cenA <- scA$centre[pop]
  sclA <- scA$scale[pop]
  
  # Compute mean SST for scenarios (use raw SST like in Code 2)
  current_sst_raw <- mean(dsub$sst_aut_raw[dsub$year %in% 2015:2020])
  future_sst_raw  <- mean(dsub$sst_aut_raw[dsub$year %in% 2035:2040])
  
  # Scale them
  current_scaled <- (current_sst_raw - cenA) / sclA
  future_scaled  <- (future_sst_raw  - cenA) / sclA
  
  # Predict
  p1 <- predict_survival(modA, vcvA, current_scaled)
  p2 <- predict_survival(modA, vcvA, future_scaled)
  
  # Save
  results_out <- rbind(
    results_out,
    data.frame(
      species = sp,
      colony  = col,
      season  = "Autumn",
      scenario = c("Current","Future"),
      surv_mu = c(p1$mu, p2$mu),
      surv_sd = c(p1$sd, p2$sd)
    )
  )
  
  # -------------------------------------------------------------------
  # WINTER model (exactly same but using sst_win)
  modW <- lm(phi_logit ~ sst_win, data = dsub, weights = phi_weight)
  vcvW <- vcov(modW)
  
  cenW <- scW$centre[pop]
  sclW <- scW$scale[pop]
  
  current_sst_win_raw <- mean(dsub$sst_win_raw[dsub$year %in% 2015:2020])
  future_sst_win_raw  <- mean(dsub$sst_win_raw[dsub$year %in% 2035:2040])
  
  current_scaled <- (current_sst_win_raw - cenW) / sclW
  future_scaled  <- (future_sst_win_raw - cenW) / sclW
  
  p1 <- predict_survival(modW, vcvW, current_scaled)
  p2 <- predict_survival(modW, vcvW, future_scaled)
  
  results_out <- rbind(
    results_out,
    data.frame(
      species = sp,
      colony  = col,
      season  = "Winter",
      scenario = c("Current","Future"),
      surv_mu = c(p1$mu, p2$mu),
      surv_sd = c(p1$sd, p2$sd)
    )
  )
}

# Final output
results_out
