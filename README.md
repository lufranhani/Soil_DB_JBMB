[![DOI](https://zenodo.org/badge/991421114.svg)](https://doi.org/10.5281/zenodo.15540526)

# Soil_DB_JBMB
Repository of results obtained by soil analysis in Bauru Botanical Garden, SP, Brazil. 
Variables considered on this study: 
Sand, Silt and Clay proportions according to Bouyoucos Method (g kg-1). 
Borum (B), Cooper (Cu), Zincum (Zn), Iron (Fe), Phosphorous (P), Sulfur (S) and Manganese (Mn) = (mg.dm-3).
Calcium (Ca), Magnesium (Mg), Potassium (K), Aluminum (Al), Potencial Acidity (H+Al), Sum of Bases (SB), and Cation Exchange Capacity (CTC) = (mmolc.dm-3). 
Nitrogen (N) and Carbon (C) = (%).
Carbon Stable Isotope (C13) and Nitrogen Stable Isotope (N15) = (δ X ‰).
Correlation of C divided by N (C/N),
Latitude (Coord_X), Longitude (Coord_Y) and Altitude (Coord_Z). 
Bulk Density (BD). 
Carbon Stocks (Cstock) = (Mg C ha -1). 

For the entire database found in this repository:
SSF = Semideciduous Seasonal Forest; 
DWS = Densely Woodded Savanna; 
DA = Deforested Area.

Sheet 1 "df" = Cstock_Total refer to the total Soil Organic Carbon Stock found as considering the entire soil profile measured (0 - 100 cm), for everyone of the 30 sampling points.

Sheet 2 "df1" = Analysis performed on the soil that correspond for the 0 - 20 cm soil layer. 
Sheet 3 "df2" = Analysis performed on the soil that correspond for the 20 - 40 cm soil layer. 
Sheet 4 "df3" = Analysis performed on the soil that correspond for the 40 - 60 cm soil layer. 
Sheet 5 "df4" = Analysis performed on the soil that correspond for the 60 - 80 cm soil layer. 
Sheet 6 "df5" = Analysis performed on the soil that correspond for the 80 - 100 cm soil layer. 
Sheet 7 "Plot" = Refer to the dataset used to create the boxplot for visualization of the results found when comparing treatments and soil depths. 

Script: "Geostatistic_Sensitivity_Analysis" used for geostatistic analysis, using "Soil_DB" database. 
Script: "Supp_Material_Papper" used to conduct pH and 13C analysis, using "Soil_DB" database. 
Script: "Pearson" used to conduct Pearson Correlation Analysis between soil properties, using an alternative Database, which are also included in the "ZIP" package. 
