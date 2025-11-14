###############################
# Ocean warming analysis of seabirds 
###############################
Corresponding author: Kate Layton-Matthews
DOI: 


We used data of sea surface temperature (extracted from seabirds' seasonal distributions), adult capture-mark-recapture and breeding success data from 26 populations of 5 seabird species (see Supplementary Information S1). 

###############################
## Data 
###############################

1) Annual survival rates for estimating and predicting non-breeding season SST effects (Data_Surv.csv)
NB: we provide  derived annual survival estimates used in all analyses,with SE. These datasets fully reproduce all results and figures in the manuscript. The raw CMR data are available upon reasonable request from the data owners at NINA. 
2) Annual breeding success values (Data_BS.csv) for estimating and predicting breeding season SST effects
3) Sea surface temperatures (SST) for the breeding season (Data_SSTbs.csv)
4) Sea surface temperatures (SST) for autumn and winter (Data_SSTaut/win.csv)

###############################
## Modelling  
###############################
1) Effects of SST on breeding success were estimated in script 1 (Model_AnnSurv_SST.R)
2) Effects of SST on survival were estimated in script 2 (Model_BS_SST.R)
3) A matrix projection model (script 3) was used to forecast population growth rates under current average SST (period 2015–2020) and future average SST (period 2035–2040) PopProj_SST.R
	--> All predictions are provided to run script three in DemRate_SST_Predictions.rdata


##############################
## Files
##############################
Main code:
1) Scripts/Model_AnnSurv_SST.R
2) Scripts/Model_BS_SST.R
3) Scripts/PopProj_SST.R

Inputs:
- Data/Data_Surv.csv
- Data/Data_BS.csv
- Data/Data_SSTbs.csv
- Data/Data_SSTaut.csv
- Data/Data_SSTwin.csv
