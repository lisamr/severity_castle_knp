# Landscape Fire Severity Analysis

## Project ideas

1.  What are the drivers of fire severity in KNP + Castle fires?
2.  Does controlling for fire direction change these drivers?
3.  How does the effect of forest structure and mortality decay with distance? How does it interact with weather (question below)?
4.  Under what weather thresholds does forest structure and mortality not matter?

## Approach

-   Gather all the relevant predictors and response variables as raster layers
-   Response:
    -   [x] fire severity. Use ML-derived metric of CBI using the Parks et al. 2019 code.
-   Possible predictors:
    -   Fire history
        -   [x] severity of last fire,
        -   [x] years since last fire.
        -   For severity---take the FRAP polygons, flatten fires, group by year into multipolygons. You'll have \~35 features, 1 for each year from 1984-2021. Use that for the GEE code (Parks et al. 2019), which won't have to be modified.
    -   topography
        -   [x] elevation
        -   [x] TPI
        -   [x] HLI
        -   [x] slope
        -   [x] aspect
    -   daily fire progression perimeters
        -   [x] from IR flights
        -   [ ] MODIS? to fill in missing days?
    -   [ ] Fire direction---heading, backing, flanking informed by aspect, slope, direction of nearest fire progression perimeter
    -   [ ] Forest structure--density, CHM where available?
    -   forest mortality
        -   [ ] Yan's model--based on 2020 NAIP
        -   [ ] eDart in the meantime--cumulative mortality from 2015-2020?
    -   weather
        -   [ ] RH
        -   [ ] wind
        -   [ ] inversions?
-   Stack these rasters---same crs, resolution, extent
-   Model the data
    -   Bayesian regression model using a distance-weighted coefficient for landscape mortality. Use \~25% of data, randomly selected, account for SA with NNGP.
    -   Run a ML model using the same predictors as above, except there's no functional structure to it.

## Repo notes

I didn't upload most of the spatial datasets bc of space. check out gitignore for that. 

## Methods notes

### Fire severity data

I used Parks et al (2019) to derive CBI from landsat imagery. That's gonna be my response variable. I also used that code to get the severity of the last fire under the KNP and Castle fire footprints.

### IR data

I downloaded IR data from <https://ftp.wildfire.gov/public/incident_specific_data/calif_s/> for the KNP and SQF fires. SQF included Castle, Shotgun and halfway through, the Rattlesnake fire. I clipped it to just include Castle fire since that's the only one I have complete data for.

I downloaded the data directly from the website using the package RCurl. The data wasn't clean and it varied whether it was exported as a shapefile, kmz, or gdb file. The shapefiles for both fires had the fewest mistakes. General steps were the following:

Download data for each day. Assume flights were taken around midnight of that date.

-   Turn into an sf object
-   Bind rows for each fire
-   Group by day (derived from createDate or filename), union so you have one multipolygon per day
-   Clip fires to castle and knp perimeters from FRAP
-   I used the FRAP perimeters for getting fire severity from GEE
-   If there were any geometrycollections, pull out just the polygons. GCs don't play nice with geometry functions.
-   Ensure that the previous day of fire was nested within that day's fire. Did this by unioning the geometry of polygon[t] with polygon[t-1] in a for-loop. Had to remove any linestrings that were created in the loop, keeping only polygons.
-   Save both versions of the data: IR polygons with and without unioning (e.g. castle_byDay.geojson, castle_byDay_u.geojson)

## Things to read

[Fitting scale-dependent landscape effect with greta - biologyforfun (lionel68.github.io)](https://lionel68.github.io/biological%20stuff/r%20and%20stat/fitting-scale-dependent-greta/)
