# Tetrahymena
Video tracking and analysis of movement of Tetrahymena pyriformis

This repository contains the code (mostly R) associated with the manuscript
https://doi.org/10.1101/2023.07.25.550219 

The data are also in the repository, with the exception of a large dataset available on figshare (10.6084/m9.figshare.23734089).

The video_tracking folder contains the python scripts used for tracking Tetrahymena cells from videos recorded under the microscope. If the scripts are used on different videos it might be necessary to change file paths and also some parameters such as the frame rate, and the scale for convering video pixels to units of length or volume (which depend on the resolution of the microscope). There is a readme file inside the folder.

The output of the video tracking are files with information about each "particle" over different frames. The header line and first data line of one such file are shown below:

"fileName", "particle", "medianFrame", "medianArea", "medianSpeed", "trajectoryLength", "autocorrelationTime", "autocorrelationTimeFrame", "autocorrelationTimeFrame2", "autocorrelationValue", "meanAbsAngle", "minSpeed", "maxSpeed", "percentile75Speed", "percentile25Speed", "meanSpeed", "meanArea", "minArea", "maxArea", "percentile75Area", "percentile25Area", "minAEllipse", "maxAEllipse", "meanAEllipse", "medianAEllipse", "minBEllipse", "maxBEllipse", "meanBEllipse", "medianBEllipse", "medianElongation", "nFrames", "meanNParticlesPerFrame", "meanDensity"

"./initial_speed_responses_2021_02_22/a_20_100_1_tested_10C_diluted10cells_ml_rep1_trajectories.csv", 0, 112.5, 592.732103692195, 345.050153900855, 224, 0.4, 12, 11.5, 0.996707680921053, 0.171420467710599, 265.154203612464, 408.394527122387, 358.185695053792, 326.71562169791, 343.508876580095, 586.683816919826, 341.269999095506, 769.781952847007, 631.221201334545, 549.111126364198, 17.6834410893635, 26.1038036987876, 21.1476598685729, 21.0407123799652, 25.6575267458534, 58.5347725221707, 41.1238528489312, 40.9340172641884, 1.94505349660619, 849, 1.26383981154299, 1160.71296255088

It is possible to get information about how each measurement is obtained by looking into the code in the video tracking folder itself, but essentially there are columns that identify the object being tracked ("filename", "particle", "trajectoryLength"), columns that describe the characteristics of the culture or of the trajectory ("meanNParticlesPerFrame", "trajectoryLength"), columns that describe the shape of the tracked object ("medianArea", "maxArea", "meanAEllipse", "meanBEllipse", "medianElongation", etc.) and columns that describe the movement ("medianSpeed", "autocorrelationTime", "meanAbsAngle").

The main data files (all with extension .csv) are:
1) /acute_speed_response/track_analysis_individual_particles.csv # a file with speed measurements of particles measured in an acute response to various temperature conditions.
2) /longer_term_speed_response/track_analysis_individual_particles.csv # same as above, but for longer term exposure to the new temperature condition

3) /Metabolic_rate_respiration/post_data/* # various files with respiration data measured in the Sensor Dish reader at different temperature conditions. Files having the extension .xls are time series of oxygen, while csv files contain information about the start and end of the recording
The file summary_MR_post_adaptation.csv contains measurements of respiration in the following form:
"intercept","slope","r2","p.val","code","temperature","concentration","treatment","replicate","density","rate.per.cell","rate.per.cell.nW"
-175.265777137724,-1.43794034568846,0.998278013581225,0,"100_15_A_8000",12.5,100,"15","A",8000,-1.79742543211058e-07,1.43367545456594

Where the most relevant columns are those that define the experiment ("temperature", "concentration", "treatment", and "replicate") and the measured rates of respiration ("rate.per.cell" and "rate.per.cell.nW"). Details of how the measurements are obtained can be found by looking into the code 

There are also other data files:
4) /population_growth_during_adaptation/growth_data_during_adaptation.csv # which contains the measurements of population density at each count and for each culture. The exact culture condition can be found in the "name" column (e.g. 20_100_1 means the culture kept at 20 degrees, 100% growth medium, replica 1). 

5) /preliminary_data_nutrients_and_population_growth/data_growth_Rate_vs_medium_concentration.csv in which cells were moved directly from 100% growth medium and 20 degrees to a range of conditions (with different growth medium concentrations and different temperatures) and their growth rate was measured (column "Growth_rate")
6) /preliminary_data_oxygen_concentration/data_cell_volume_vs_oxygen.csv # where there are measurements of cell volume for cells kept at different oxygen concentrations
7) /speed_over_time_preliminary_experiment/track_analysis_results_individual_particles_preliminary_experiment_long_tracking.csv # with measurements of acute speed responses from preliminary experiments and over a relatively long period of time (20 minutes)

The R code in each folder is meant to run on the data file in the same folder (or sometimes on data from different folders, but the paths indicate which files are read by each script) and then each script produces a number of figures or processed additional datasets (e.g. tables), which are saved to disk.


