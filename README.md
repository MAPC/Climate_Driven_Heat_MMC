# Climate_Driven_Heat_MMC
Repository storing tools for Metro Mayors Climate-Driven Heat analysis and data development.

Task 2.1: Development of a Massachusetts-wide Land Surface Temperature (LST) index raster dataset.
The LST index is derived from Landsat 8 surface temperature estimates. Status: Beta testing.

Task 2.2: Compares spatial patterns in air temperature within Metro Mayors Coalition (MMC) communities, 
as measured by daily maximum air temperature at Weather Underground Personal Weather Stations (PWS) 
within MMC communities, with the spatial patterns of LST index at each PWS. Contains tools to extract
land cover composition and average LST index within a range of buffer diameters of the PWS locations,
and documents linear model selection and calibration process producing a simple linear regression model
relating the difference between spatial pattern of air temperature vs. LST index with local land cover
composition. Produces a vector-format polygon layer which documents where in the MMC region air temperatures
can be expected to be higher, lower, or about the same as expected based on the spatial pattern of 
LST index. Status: Draft.

Task 2.3: Contains tools for modelling the relationship between regional average air temperature and 
localized air temperature as measured by Weather Underground PWS with land cover composition metrics
as key covariates. Produces spatially explicit summaries of extreme heat intensity-duration-frequency 
metrics by combining regional temperature projections with aforementioned regional air temperature/local 
air temperature statistical model and future land cover scenarios. Status: Under development.

Task 2.4: Updates the MAPC Climate Vulnerability Index (https://climate-vulnerability.mapc.org/) to adopt 
the LST index data produced in Task 2.1 as the heat exposure indicator, replacing the 2016 single-date land 
surface temperature raster dataset which is the current extreme heat exposure indicator. Status: Under 
development. 


