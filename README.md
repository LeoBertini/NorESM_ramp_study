## This repository contains source code used in the following study: 

**Biogeochemical Timescales of Climate Change Onset and Recovery in the North Atlantic Interior Under Rapid Atmospheric CO2 Forcing**

Leonardo Bertini, Jerry Tjiputra
First published: 23 March 2022

Journal of Geophysical Research: Oceans
https://doi.org/10.1029/2021JC017929

## A description of each script is given below:
**Disclaimer**: CODE IS PUBLISHED AS IS. 

### 1 - Vertical interpolation
[pre_ind_leo_grunch_BIOGEOCHEM.m](https://github.com/LeoBertini/NorESM_ramp_study/blob/main/pre_ind_leo_grunch_BIOGEOCHEM.m')  
Vertical interpolation code for transforming isopycnal fields into depth fields. 

### 2 - Annual means (Biogeochemical)
[PI_BGC_trends_annual_averages.m](https://github.com/LeoBertini/NorESM_ramp_study/blob/main/PI_BGC_trends_annual_averages.m')  
Routine to obtain the annual average for Biogeochemical fields (Dissolved Oxygen, pH, Omega calcite, Export Production).

### 3 - Annual means (Physics)
[PI_annual_averages_physics.m](https://github.com/LeoBertini/NorESM_ramp_study/blob/main/PI_annual_averages_physics.m')  
Routine to obtain the annual average for Physical fields (Temperature and Salinity).

### 4 - Apparent Oxygen Utilization
[AOUcalculation.m](https://github.com/LeoBertini/NorESM_ramp_study/blob/main/AOUcalculation.m)  
Calculation of Apparent Oxygen Utilization (AOU) 

### 5 - Timescales 
[Time_scale_analyses_BIOGEOCHEM.m](https://github.com/LeoBertini/NorESM_ramp_study/blob/main/Time_scale_analyses_BIOGEOCHEM.m)  
Timescale analyses of Tracertmax, Time of Departure (ToD) and Time of Recovery (Trec)

### 6 - Stats
[Stats_compare_PI_vs_Ramps.m](https://github.com/LeoBertini/NorESM_ramp_study/blob/main/Stats_compare_PI_vs_Ramps.m)  
Statistical analyses comparing the last 30 years of the Pre-Industrial simulation against 30-year windows: 
 - 1) end of the Mitigation phase (years 250-280), 
 - 2) the middle of the Extension phase (years 350-380) and 
 - 3) end of the Extension phase (years 450-480) 

[Code for Figures](https://github.com/LeoBertini/NorESM_ramp_study/tree/main/Code_for_figures)  

- Vertical profiles of Timescales across the entire NAtl (0-65N)
- Time-evolution of biogeochemical drivers at depth for the entire NAtl per lat band
- Animations


## DESCRIPTION OF POST-PROCESSED DATA

### [Results: Statistical tests](https://drive.google.com/drive/folders/1o62Zr4Tw-npM6_kNsQUzc2U5cDyGH1_T?usp=sharing)  

Contains percentage change and the statistical comparisson  results between Pre-Industrial and: 
1) years 250-280; 
2) years 350-380 and 
3) years 350-380

Files are for each variable : pH, AOU, o2sat (oxygen saturation), o2 (dissolved oxygen), omegac( Calcite Saturation State), templvl (temperature), epc100 (POCexport)


### [Results : Time scale analyses](https://drive.google.com/drive/folders/1M0vbH0hiiDIqOshxMqmSYp9oZ_Sb1Sd9?usp=sharing)  

Contains Time scale estimates for Time of Departure (ToD), Time of Recovery (Trec) and extent of change (Maxchange) 
for different prescribed envelopes of 1 standard deviation, 2 standard diviations and 3 standard deviations. 

Obs. Trec estimates are based on a linear regression using the last 100 years of the time series whenever there is no recovery whitin the time span of the simulation.  



#__Contact :__
Leonardo Bertini
leonardo.bertini@imbrsea.eu | l.bertini@nhm.ac.uk 
