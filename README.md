# Mburyatense_data
This project contains the raw data files and MATLAB scripts needed to read these files for the paper "A resource allocation model for the Methanotroph Methylotuvimicrobium buryatense".

The data of three experimental runs has been used in the paper. The "read_data" MATLAB scripts are used to read and plot the raw data of one specific experiment.

Experiment of 15/06/2020
  - Read_data_15062020.m
  - BiomassHPLC_20200615.xlsx
  - BlueVis_Export_beginning_from_2020_6_16_10_46.csv
  - CTPC06603.20200616 5GB1.Control.csv

Experiment of 14/07/2020
  - Read_data_13072020.m
  - BiomassHPLC_20200714.xlsx
  - BlueVis_Export_beginning_from_2020_7_14_10_0.csv
  - CTPC06603.20200714 - 5GB1 - KDB.Control.csv

Experiment of 27/10/2020
  - Read_Data_27102020.m
  - BiomassHPLC_20201027.xlsx
  - BlueVis_Export_beginning_from_2020_10_27_11_5.csv
  - CTPC06603.Manager 88.Control.csv

Additionally, the files Enz_list.xlsx and Model_Mburyatense3.xls contain a full list of enzymes taken into account in the deFBA model for M. buryatense and a detailed description of all included reactions and metabolites, repsectively.

The folder deFBA_Mbur contains the files needed to run the deFBA model for the M. buryatense model, the corresponding paramter estimation, the RBA analysis and the sensitivity analysis. Only the deFBA core code is missing due to copyright reasons. This code is available upon request by contacting Prof. Steffen Waldherr.

All files in this repository are covered by the GNU GENERAL PUBLIC LICENSE.
