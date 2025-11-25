library("popbio")
library("arm")
library("boot")
library("ggplot2")
library("pracma")
library("plyr");library("dplyr")
library("ggpattern")
library("tidyr")
library("gridExtra")
library("patchwork")
library("grid")
library("gtable")

########## Description ########## 
#' Script to calculate impact weights for climate change stressor
#' Current impact weight = historic lambda / current lambda
#' Future impact weight = historic lambda / future lambda 
#' 
#' Option to loop through each populations and saves data as a dataframe 
#' Matrices are made for N simulations of vital rates and mean (SD) is calculated 

############# Load data ############# 
SSTeff = readRDS(file="DemRate_SST_Predictions.RDS")
SPMeta = read.csv(file="LookUp_Species.csv",sep=";") # for Age at Maturity
SPMeta$Species=tolower(SPMeta$Species)

# Unique.id for all dfs
unique.id=unique(paste(SSTeff$Species,SSTeff$Colony,sep=" "))

# Pars for for loops
npop = length(unique(unique.id)) # j
nsim = 10000 # i

############# Empty lists to hold results ############# 
lambdas = data.frame(Species=NA,
                     Colony=NA, 
                     lambda_C=rep(NA,nsim),
                     lambda_F=rep(NA,nsim))

lambdas.list = lapply(1:npop, function(x) lambdas)
rm(lambdas)

############# Make matrices of different size ############# 
atpu <- matrix(0, nrow = SPMeta[SPMeta$Species=="atpu","Age.maturity"],  ncol = SPMeta[SPMeta$Species=="atpu","Age.maturity"])
blki <- matrix(0, nrow = SPMeta[SPMeta$Species=="blki","Age.maturity"],  ncol = SPMeta[SPMeta$Species=="blki","Age.maturity"])
cogu <- matrix(0, nrow = SPMeta[SPMeta$Species=="cogu","Age.maturity"],  ncol = SPMeta[SPMeta$Species=="cogu","Age.maturity"])
brgu <- matrix(0, nrow = SPMeta[SPMeta$Species=="brgu","Age.maturity"],  ncol = SPMeta[SPMeta$Species=="brgu","Age.maturity"])
liau <- matrix(0, nrow = SPMeta[SPMeta$Species=="liau","Age.maturity"],  ncol = SPMeta[SPMeta$Species=="liau","Age.maturity"])

SpMat <- list(atpu = atpu, blki = blki, cogu = cogu, brgu = brgu, liau = liau)

for (j in 1:npop){
  
  lambdas.list[[j]]$Species = substr(unique.id[j],1,4)
  lambdas.list[[j]]$Colony = substr(unique.id[j],6,9)
  
  Species=substr(unique.id[j],1,4)
  Colony=substr(unique.id[j],6,9)
  
  ## Subset data for focal population 
  SSTeff_sub=SSTeff[SSTeff$Species==Species & SSTeff$Colony==Colony,]
  SPMeta_sub=SPMeta[SPMeta$Species==Species,]
  
  ########## Demographic rate dataframe  ########## 
  
  vr_df=data.frame(Species=Species,
                   Colony=Colony,
                   Scenario=c("Current","Future"),
                   adsurv_mean=NA,
                   adsurv_sd=NA,
                   prod_mean=NA,
                   prod_sd=NA,
                   jsurv_mean=NA,
                   jsurv_sd=NA,
                   imsurv_mean=NA,
                   imsurv_sd=NA)
  
  #Collate data depending on what level is available(Population-Scenario, Species-Scenario or Species only)
  vr_df$adsurv_mean=SSTeff_sub$mean_surv
  vr_df$adsurv_sd=SSTeff_sub$se_surv
  vr_df$prod_mean=SSTeff_sub$mean_prod
  vr_df$prod_sd=SSTeff_sub$se_prod
  vr_df$jsurv_mean=SSTeff_sub$mean_surv*0.75
  vr_df$jsurv_sd=SSTeff_sub$se_surv
  vr_df$imsurv_mean=SSTeff_sub$mean_surv*0.875
  vr_df$imsurv_sd=SSTeff_sub$se_surv
  
  ############# Vital rate distributions  ############# 
  vr_C=vr_df[vr_df$Scenario =="Current",]
  vr_F=vr_df[vr_df$Scenario =="Future",]
  
  vr_C_sim = vr_F_sim=matrix(data=NA,nrow=nsim,ncol=4)
  colnames(vr_C_sim) = colnames(vr_F_sim)=c("sj","sim","sad","prod")
  
  # Simulate for Current
  vr_C_sim[,"sad"]=invlogit(rnorm(nsim,mean=logit(vr_C$adsurv_mean),sd = vr_C$adsurv_sd))
  vr_C_sim[,"prod"]=invlogit(rnorm(nsim,mean=logit(vr_C$prod_mean),sd = vr_C$prod_sd))
  vr_C_sim[,"sim"]=invlogit(rnorm(nsim,mean=logit(vr_C$imsurv_mean),sd = vr_C$imsurv_sd))
  vr_C_sim[,"sj"]=invlogit(rnorm(nsim,mean=logit(vr_C$jsurv_mean),sd = vr_C$jsurv_sd))
  
  # Simulate for Future
  vr_F_sim[,"sad"]=invlogit(rnorm(nsim,mean=logit(vr_F$adsurv_mean),sd = vr_F$adsurv_sd))
  vr_F_sim[,"prod"]=invlogit(rnorm(nsim,mean=logit(vr_F$prod_mean),sd = vr_F$prod_sd))
  vr_F_sim[,"sim"]=invlogit(rnorm(nsim,mean=logit(vr_F$imsurv_mean),sd = vr_F$imsurv_sd))
  vr_F_sim[,"sj"]=invlogit(rnorm(nsim,mean=logit(vr_F$jsurv_mean),sd = vr_F$jsurv_sd))
  
  ############# Start of vital rate sims for loop  ############# 
  
  for (i in 1:nsim){
    ########################
    # Current matrix
    ########################
    
    # Save vital rates for historic matrix 
    sj=vr_C_sim[i,"sj"]
    sim=vr_C_sim[i,"sim"]
    sad=vr_C_sim[i,"sad"]
    prod=vr_C_sim[i,"prod"]
    
    # Construct current matrix
    Amat= SpMat[[Species]]
    
    Amat[row(Amat) == col(Amat) + 1] <- sim
    Amat[2,1]=sj
    Amat[nrow(Amat),ncol(Amat)] <- sad
    Amat[1,ncol(Amat)] <- sad*prod*0.5
    
    # Save current lambda
    lambdas.list[[j]]$lambda_C[i] = eigen.analysis(Amat)$lambda
    
    ########################
    # Future matrix
    ########################
    
    # Save vital rates for historic matrix 
    sj=vr_C_sim[i,"sj"]
    sim=vr_F_sim[i,"sim"]
    sad=vr_F_sim[i,"sad"]
    prod=vr_F_sim[i,"prod"]
    
    # Construct current matrix
    Amat= SpMat[[Species]]
    
    Amat[row(Amat) == col(Amat) + 1] <- sim
    Amat[2,1]=sj
    Amat[nrow(Amat),ncol(Amat)] <- sad
    Amat[1,ncol(Amat)] <- sad*prod*0.5
    
    # Save current lambda
    lambdas.list[[j]]$lambda_F[i] = eigen.analysis(Amat)$lambda
    
  } # i end vital rate simulation loop
} # j end population loop 

# ~~~~~~~~~~~~~~~~~~~~~~ # END for loop for vital rates # ~~~~~~~~~~~~~~~~~~~~~~ # 

###############################
# Collapse lists into single dataframe and calculate summary statistics 
##############################

lambdas.df=do.call(rbind,lambdas.list)

########## Calculate mean lambda per population and scenario (mean, SE & 95% CIs) 
lambdas.df.long=pivot_longer(lambdas.df,cols = c( "lambda_C", "lambda_F"),
                             names_to = "Scenario",
                             values_to = "Lambda")
lambdas.df.long$Scenario=factor(lambdas.df.long$Scenario,levels=c("lambda_C","lambda_F"))
levels(lambdas.df.long$Scenario) = c("Current","Future")

lambdas.df.long.sum=lambdas.df.long %>%
  group_by(Species,Colony,Scenario) %>%
  summarise(Mean = mean(Lambda),
            sd = sd(Lambda),
            lcl = quantile(Lambda, probs = 0.20),
            ucl = quantile(Lambda, probs = 0.80))

print(lambdas.df.long.sum,n=nrow(lambdas.df.long.sum))
