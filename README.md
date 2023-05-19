# Landscape Fire Severity Analysis

## Project ideas

1.  What are the drivers of fire severity in KNP + Castle fires?
2.  Does controlling for fire direction change these drivers?
3.  How does the effect of forest structure and mortality decay with distance? How does it interact with weather (question below)?
4.  Under what weather thresholds does forest structure and mortality not matter?

Reasons to look at scale of mortality: managers need to know where to conduct fuel treatments. Say they are thinning in sequoia groves, they will need to thin around the grove too to make it effective. How big of a buffer do they need?

## Approach

-   Gather all the relevant predictors and response variables as raster layers
-   Response:
    -   [x] fire severity. Use ML-derived metric of CBI using the Parks et al. 2019 code.
        -   [x] figure out why CBI is distorted on the upper range. \~3% \>2.25 in KNP fire, maybe a little higher for Castle. I didn't include Landsat9 to estimate KNP, but that wouldn't address why Castle fire is also low. Using values from MTBS for now
        -   [x] differences between deprecated and current Landsat datasets needed scale and offset corrected. GEE code fixed now.
-   Possible predictors:
    -   Fire history
        -   [x] severity of last fire,
        -   [x] years since last fire.
        -   For severity---take the FRAP polygons, flatten fires, group by year into multipolygons. You'll have \~35 features, 1 for each year from 1984-2021. Use that for the GEE code (Parks et al. 2019), which won't have to be modified.
    -   topography
        -   [x] elevation
        -   [x] TPI -- need to buffer dem by 1000m since your window goes up to that.
        -   [x] HLI
        -   [x] slope
        -   [x] aspect
        -   [x] add hillshade for viz
    -   daily fire progression perimeters
        -   [x] from IR flights
        -   [x] VIIRS to fill in missing days? use it completely?
        -   [ ] for now use viirs. but manually assess when IR needs to be replaced with IR.
    -   [x] Fire direction---heading, backing, flanking informed by aspect, slope, direction of nearest fire progression perimeter
    -   [x] Forest structure--density, CHM where available? SD of NDVI might work eg. Koontz EcoLett 2020?.
        -   [x] added NDVI, should be good enough. used composited landsat imagery year of fire
            -   [x] make a wider buffer around fire when exporting. 5km?
    -   forest mortality
        -   [ ] Yan's model--based on 2020 NAIP
        -   [ ] ~~eDart~~ F3 in the meantime. Looked at F3 mortality dataset. it's from MMI. need to compile that from individual years.
    -   weather
        -   [x] RH--gridmet?
        -   [x] ~~wind~~
        -   [ ] inversions?
-   [x] Stack these rasters---same crs, resolution, extent
-   Model the data
    -   Bayesian regression model using either a ~~distance-weighted coefficient~~ scale-selecting coef for landscape mortality. Use \~5%? of data, randomly selected, account for SA with NNGP.
        -   [NNGP in Stan](https://mc-stan.org/users/documentation/case-studies/nngp.html). I've tested this out before with good success.

## Repo notes

-   I didn't upload most of the spatial datasets bc of space. check out gitignore for that.

-   at this point, still compiling much of the data. Waiting on NAIP-derived mortality layer from collaborator, need to get daily fire progression data. IR flight data has missing days. Might use VIIRS satellite data to interpolate boundaries. See: <https://www.mdpi.com/2072-4292/12/12/2061>

-   I'm testing out how to estimate the scale at which mortality is important. Currently using sd of NDVI to test Stan code. There are 2 approaches--1) scale-selecting (SS) approach determines distance (with uncertainty) at which effect is greatest [(scripts/03_statistical_models/xx_test_ssmodel.R)](https://github.com/lisamr/severity_castle_knp/blob/main/scripts/03_statistical_models/xx_test_ssmodel.R) and 2) distance-decay (DD) approach assumes that effect decays with distance, determined by a kernel function, and estimates parameter(s) controlling that decay [(scripts/03_statistical_models/xx_test_ddmodel.R)](https://github.com/lisamr/severity_castle_knp/blob/main/scripts/03_statistical_models/xx_test_ddmodel.R).

    -   SS summarizes data in buffers and DD summarizes data into non-overlapping rings. More info on these approaches in the linked literature down below.

    -   I've found a way to do both approaches in Stan and both work well now. The D-D approach seems to require a bit more data, but will probalby be fine.

## Methods notes

### Fire severity data

I used Parks et al (2019) to derive CBI from landsat imagery. That's gonna be my response variable. I also used that code to get the severity of the last fire under the KNP and Castle fire footprints.

-   [x] About that: the default dates are 6/1 to 9/30. creek fire started 8/26. double check that the pre-fire imagery is for the previous year, not the same year as fire.
-   [x] Parks code uses deprecated Landsat datasets (collection1 instead of current collection2). Using the current Collection 2 dataset, which includes Landsat 9, you need to correct code for differences in scale and offset. GEE code corrected and datasets in this repo are updated.

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

-   [Fitting scale-dependent landscape effect with greta - biologyforfun (lionel68.github.io)](https://lionel68.github.io/biological%20stuff/r%20and%20stat/fitting-scale-dependent-greta/)
    -   The code in that blog post isn't that useful, but it was an entry way to the methods from Miguet et al 2017.
-   [Miguet 2017 paper](https://doi.org/10.1111/2041-210X.12830) explaining how to estimate the decay of distance-dependent effects
-   [Moll et al. 2020](https://doi.org/10.1111/ecog.04762) on how urbanization impacts wildlife An application of the Miguet paper, including JAGS code
    -   I translated these models into Stan. I tested it with simulated data (assumed that rings of SD of NDVI predicted normally distr. response variable) and it worked well. When I tried on real data, the scale parameter wasn't informed by the data whatsoever. Maybe I need more data?
    -   I just needed more data. Everything is working well now.
-   [Monroe et al. 2022](https://doi.org/10.1002/ecs2.4320) used a scale-selecting parameter with linear interpolation to estimate the scale of effect (and uncertainty around it). Paper was on sage grouse counts and surrounding sagebrush habitat. I chatted with lead author a bit over Teams and he's readily communicative which is a plus
    -   Pretend you have a whole bunch of focal sites and you created multiple concentric buffers around each one. These buffers extend from 50 to 2000 in 50m intervals and will be used to estimate which scale of the effect is strongest. Now you have a data matrix $\mathbf{E}$ containing the mean tree mortality within each of these buffers with N rows for each focal site and M columns for each scale. In the model, you estimate a parameter $\varsigma$, bounded between 1 and M, to select which scale of $\mathbf{E}$ is strongest in predicting your response variable. Because the data is computed in discrete steps, you use linear interpolation to estimate the effect of $\mathbf{E}$. E.g. if $\varsigma$ is 3.4, its calculated as $\mathbf{E[,3]} \times .6 + \mathbf{E[,4]} \times .4$.
    -   Internally, to do this interpolation, you need to calculate the floor of $\varsigma$ and treat it as an integer in order to index $\mathbf{E}$. Stan doesn't allow parameters to be integers. I found this workaround on the Stan forum of someone else asking the exact question I needed [here](https://discourse.mc-stan.org/t/workaround-using-truncated-parameters-for-indexing/28971/2).
    -   I tested out this approach and so far it's working, even with real data.
