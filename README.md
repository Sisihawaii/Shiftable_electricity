# Shiftable_electricity
Code for paper: Using Temperature Sensitivity to Estimate Shiftable Electricity Demand
Co-authors: Michael Roberts, Eleanor Yuan, James Jones, Matthias Fripp

0: Preparing the data (from James Jones previous work):
0.1. Download data (0.1. Download.R)
  1A: Download NARR
  1B: Download EIA 930
  1C: Download BA shapefile
  1D: Download population grid
  1E: Download EIA Meta - BA name match file
  
0.2. Clean Data
  2A: Reproject BA shapefile to LCC to match NARR (0.2. reprojSF.R)
  2B: Match EIA ids to BASF ids (0.2. matchIDs.R)
  2C: Consolidate BAs in each interconnection (0.2. mergeBAs.R)
  
0.3. Calculate Weather
  3A: Calculate pop weighted weather for each consolidated BA (0.3. calcBAWDD.R)
  
--------------------------------------

Main code:

1. climavg_kaholo.R
  - Calculate 12 years average climate for cdh and hdh per grid cell
  - Interact climate averages with cdh and hdh per cell respectively
  - Pop weighted for each BA

2. data_cleaning.R
  - Use clean data from https://github.com/truggles/EIA_Cleaned_Hourly_Electricity_Demand_Code
  - Adjust clean data to fit the regression and combine with weather data

3. trainBAs.R
  - Train the optimal splines using CV for each BA in East and West IC

4. reg_e15ba.R
4. reg_w15ba.R
  - Run regression to calculate flexible load with the optimal splines for east and west IC

5. Regcoef_visual.R
  - Visualization on range of cdh/hdh coef for each BA for different models
  - RMSE reduction between different models

6a. flattendata.R
  - Flatten demand across the BA regions' individual timezone if consecutive 23 or 24 hours are present 
  - Find residual between estimated and real demand data.

6b. summarytable.R
  - Summarise daily or annual % peak reduction, % base increase, % sd/cv reduction and daily % count of completely flattenable days (sd = 0).
  - Values calculated at different levels of flexibility (alpha = 0, 0.25, 0.5, 1), different levels of grid transmission (within each BA region, each interconnect, the nation), and across all seasons generally and for Winter and Summer months specifically.

6c. alphashifts.R
  - Generate values and graphs to showcase total % sd reduction for different levels of flexibility (alpha = 0:1) for each BA region daily and annually and for each season.

6d. SDred_map.R
  - Visualization map on daily demand SD reduction for the US

7. Climate2C_kaholo.R
  - Calculate new climate data with +2C temperature increase everywhere
  - Run the regression model with new climate data and predict the electricity demand/flexible load under climate change

8. remodel_2c.R
  - Rerun the model with new climate data without climate interactions

9. Ercot.R
  - Repeat above for Ercot

10a. flattendata2c.R
  - Same as flattendata.R but applied for +2C temperature increase data and including residual values calculated from flattendata.R.

10b. basepeak_bargraph.R
  - Generate graph to showcase daily or annual peak and base demand in relation to average demand of the BA region across different levels of flexibility and different levels of grid transmission
  - Whiskers showcasing 1st and 99th percentile are also plotted.

11. appendix_seascdhhdh_temp.R
  - Generate graphs to showcase average demand for each season in each BA region on a Demand Index scale (demand divided by average demand of the BA region) according to the temperature and a bar graph below each line plot to showcase total count of hours that recorded each temperature in that BA region.
  - Generate graphs to showcase average demand index value according to if the hour is counted as CDH or HDH and a bar graph below each line plot to showcase total count of hours recorded for each CDH or HDH value (above or below 18C), bars are overlayed.

--------------------------------------

Function code:

- function_ba.R: functions for 0. data preparation 
- functions.R: all the functions used in this work
- spec_chart_function.R: 
  # Authors: Ariel Ortiz-Bobea (ao332@cornell.edu).
  # Version: March 10, 2020


