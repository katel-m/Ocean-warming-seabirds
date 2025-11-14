###############################
# Ocean warming analysis of seabirds 
###############################
Corresponding author: Kate Layton-Matthews
DOI: 

## Species

We used data of sea surface temperature (extracted from seabirds' seasonal distributions), adult capture-mark-recapture and breeding success data from 26 populations of 5 seabird species (see Supplementary Information S1). 

###############################
## Data provided are 
###############################

1) Annual survival rates for estimating and predicting non-breeding season SST effects (Data_Surv.csv)
2) Annual breeding success values (Data_BS.csv) for estimating and predicting breeding season SST effects
3) Sea surface temperatures (SST) for the breeding season (Data_SSTbs.csv)
4) Sea surface temperatures (SST) for autumn and winter (Data_SSTaut/win.csv)

## Model
Effects of SST on survival and breeding success were estimated used scripts 1 and 2. 
A matrix projection model (script 3) was used to forecast population growth rates under current average SST (period 2015–2020) and future average SST (period 2035–2040).


##############################
## Files Overview
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

Outputs: 
Predictions are demographic rates saved in DemRate_SST_Predictions.rdata
