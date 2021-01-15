# MATLABGravityandWaterLevels
This code was primary created for the analysis of changing head values to determine likelihoods.
The change in Water Levels and the gravity of the cells is calculate in this code. Run the water_level_means_other.m
first for the desired ensemble set (See Modelgenerationcode repository to find how the values are created). After
the likelihoods and POI values created by the first code can be used in UtilityModel_May to create utilitys of each model
and weight the values to find the errors and the value of the set observations. The base_top_elev.csv can be altered but it 
is included to show what was used in the initial analysis.