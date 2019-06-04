# metHIEV
process met data from EucFACE and upload to hiev. 
The code takes met data (CO2, PAR, Tair, PPT, windspeed,and RH) from EUCFACE and aggragate to 30 min timestep. The met data are the average of all rings. Co2 are two means of elevated and ambient rings,respectively. 

The automet is the file that is used to generate data uploaded for hiev automatically.

Oldmet is the one that processes the previous met data. It uses the same functions as auto met, which are stored in the function file.
