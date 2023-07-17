# Landscape Fire Severity Analysis

## Project ideas

1.  What are the drivers of fire severity in KNP + Castle fires?
2.  How does the effect of forest structure and mortality decay with distance?
3.  Under what weather thresholds does forest structure and mortality not matter?

Reasons to look at scale of mortality: managers need to know where to conduct fuel treatments. Say they are thinning in sequoia groves, they will need to thin around the grove too to make it effective. How big of a buffer do they need?

## Approach

Using remote sensing data, we will analyze fire severity with respect to surrounding tree mortality, while controlling for relevant abiotic and biotic variables. Model development has begun on the KNP fire, but a larger goal is to expand the analysis to many fires across the Sierra Nevada (see section on Spatial Autocorrelation). A distance decay kernel will be embedded within a Bayesian model to estimate the landscape-level effect of tree mortality. We will then explore how different thinning intensities and treatment areas alter predicted fire severity with a focus on giant sequoia groves.

![A. Landscape effects of tree mortality on fire severity (cells) will analyzed by accounting for dead trees far from a focal site, while assuming that the importance (shading) decreases with distance. B) The scale that contains 90% of the mortality effect. Here, mortality is most important within 0.6 km. C) We expect fuel treatments over larger areas and of higher intensity to reduce risk of high severity fire.](figures/calfire_final.png)

Fire severity $y_i$ at observation $i$ is categorized into 4 ordered levels--unchanged, low, moderate, and high--and is modeled with an ordered logit likelihood, where $\phi_i$ controls the ordinal probabilities and $\kappa$ estimates the internal cutpoints (equivalent to the intercepts of each ordered category). For more information on ordinal regressions in Stan, see this tutorial by Michael Betancourt [here](https://betanalpha.github.io/assets/case_studies/ordinal_regression.html).

$$
y_i \sim \text{ordered-logit}(\phi_i, \kappa) \\
\phi_i = XB + (X_{DD} W)\beta + \text{interaction terms} \\
w_{0[m]} = exp(-2*distance_m^2 / \delta^2) * area_m \\
W = \frac{w_0}{\sum^M_m w_{0[m]} }
$$

The equation for $\phi$ is where you'd include your typical environmental predictors, including a matrix of variables in $X$ and vector of parameters $B$. Landscape-level effects are estimated by weighting landscape variables, which have been discretized into $M$ non-overlapping concentric rings, by a function that decreases with distance and standarized to sum to one (see equation $W$; for more information consult [Miguet et al. 2017](https://doi.org/10.1111/2041-210X.12830) and [Moll et al. 2020](https://doi.org/10.1111/ecog.04762)). Common distance weighting functions including negative exponential or negative quadratic exponential (as shown above) and can be assessed for predictive power using model comparison. The smoothing function assumes that the effect of a landscape variable is strongest close to the focal site and decays as a function of distance. In the equation above, $\delta$ controls the length scale of the decay function and coursely describes the distance (in km) at which the effect is about half. Density of dead trees and live trees are both included as distance dependent effects. Distance dependent effects are included in the model as an array of matrices with N rows and M columns, where elements are the mean of the variable within each $i$ focal site's ring $m$. The matrices $X_{DD}$ are multipled by $W$, the weighting function and the overall effects of the landscape varialbes are determined by $\beta$. Interactions between variables, including the distance dependent varialbes, can also be encoded within the equation for $\phi$. To further understand the model structure, inspect the full Stan model here: "severity_castle_knp/scripts/stan_models/ordinal_DDarray_int.stan".

## Spatial Autocorrelation

I sampled data from the KNP fire on a regular grid with a spacing of 180 m (see scripts/03_sample_data/01_assign_grid_locations.R). When testing the model (see scripts/04_statistical_models/01_ordinal_DD_mods.R) I subsampled the data so that it was sampled on a sparser grid in order to see how sampling density affected parameter estimation. I noticed that when the sampling grid was dense (e.g. 180 m), the posterior for $\delta$, the length scale parameter, was very narrow. And when I sampled on a much sparse grid, the posterior became very wide and effectively can't be estimated.

It is possible that the model fails to estimate $\delta$ as the grid becomes sparser because the model requires a larger sample size. It is also possible that when the sampling grid is tight, I'm actually detecting spatial autocorrelation (SA) instead of a distance-dependent effect.

Before going forward, it's important to test if the model can detect a distance dependent effect when SA is present and figure out at what density the sampling grid must me, while making some assumptions about the scale of SA. I would simulate data with and without SA embedded in the data generation process. Then I would run the models as is to see if if SA is causing these patterns affecting the $\delta$ posterior. That's my first gut feeling, but it will likely need further investigation.

If SA is messing things up, either SA needs to be incorporated into the model or the data needs to be sampled on a sparser grid. SA is super computationally expensive and another can of worms. I started going down that rabbit hole in the folder scripts/hilbert_space, which is a way to *as efficiently as possible* account for SA (see <https://discourse.mc-stan.org/t/practical-hilbert-space-approximate-bayesian-gaussian-processes-for-probabilistic-programming/14633>). OR you just sample on a sparser grid where SA isn't much of a concern. If you go that route however, your sample size will decrease, possibly below a critical threshold in order to adequately estimate parameters. To make up for that, I'd expand the study to multiple fires over many years (constrain to sierra nevada) and include fire ID as a varying intercept in the model. This will require more work to compile all the geospatial data, but it will also make the study more powerful as it would then be a regional/statewide study.

## Data layers

-   Gather all the relevant predictors and response variables as raster layers

-   Response:

    -   [x] fire severity. Use ML-derived metric of CBI using the Parks et al. 2019 code.
        -   [x] figure out why CBI is distorted on the upper range. \~3% \>2.25 in KNP fire, maybe a little higher for Castle. I didn't include Landsat9 to estimate KNP, but that wouldn't address why Castle fire is also low. Using values from MTBS for now
        -   [x] differences between deprecated and current Landsat datasets needed scale and offset corrected. GEE code fixed now.

-   Possible predictors:

    -   forest mortality
        -   [x] Yan's model--based on 2020 NAIP
            -   [x] check if red trees are just tracking hardwoods (e.g. buckeyes, blue oaks). if so, would want to mask those out and only consider dead conifers for both red and grey stage trees.
        -   [ ] ~~eDart~~ F3 in the meantime. Looked at F3 mortality dataset. it's from MMI. need to compile that from individual years.
        -   [x] F3 predictions of live trees. Using TPA\--live tree counts 24+inches in DBH
    -   live forest density
        -   depending on what gets used for measuring mortality, use an equivalent dataset to measure forest density (e.g. %canopy area or \# of live trees)
            -   I have live tree datasets in the data folder from Salo (canopy area) and F3 (basal area and TPA)
        -   more dead trees will exist where there are more trees. you need to control for density of live trees to isolate the direct effect of dead ones.
        -   To keep things on equal playing fields, also make live trees a distance dependent predictor. The Stan model is written to accomodate a list of distance dependent variables.
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
            -   [x] Dates: flights usually around midnight. dates recorded as next day (eg. 10/16 \@ 2300 = 10/17). adjust the date so the climate is extracted from previous day.
        -   [x] VIIRS to fill in missing days? use it completely?
            -   [x] the satellites pass over in the afternoon (12-3pm) and evening (1-3am). I adjusted the dates so that observations at night refer to the day before since most of the growth was happening when it was light out.

            -   [x] used DBSCAN to cluster points using 1125 m and 3 points min. values come from Briones 2020 paper. Then used concave hull to make perimeters. should approximate that ArcPro tool.
        -   [ ] for now use viirs. but manually assess when IR needs to be replaced with IR.
    -   [x] Fire direction---heading, backing, flanking informed by aspect, slope, direction of nearest fire progression perimeter
        -   [ ]
    -   [x] Forest structure--density, CHM where available? SD of NDVI might work eg. Koontz EcoLett 2020?.
        -   [x] added NDVI, should be good enough. used composited landsat imagery year of fire
            -   [x] make a wider buffer around fire when exporting. 5km?
    -   weather
        -   [x] might want daytime info. viirs comes by at \~noon (Day) and \~midnight (Night). Night show way more growth bc it's capturing afternoon fires.
        -   [x] climate variables from Park Williams. Metadata/explanations are in that folder (wiliams_weather)
            -   [x] daily weather when fire passed over. outputted as sf points. sampled on a 180m grid.

                -   [ ] should the weather also be resampled at a 30m grid? probably.

            -   [x] climate averages 2000-2020 May1-Oct31. raster should be resampled at a 30m grid.

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

-   [x] About that: the default dates are 6/1 to 9/30. creek fire started 8/26. double check that the pre-fire imagery is for the previous year, not the same year as fire. edit: yep.
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

-   <https://link.springer.com/article/10.1007/s10980-022-01437-5>
    -   a package for dealing with distance dependant varialbes and spatial autocorrelation. likely too much of a black box for developing an interesting model but there might be good text in there to learn from.
-   [Fitting scale-dependent landscape effect with greta - biologyforfun (lionel68.github.io)](https://lionel68.github.io/biological%20stuff/r%20and%20stat/fitting-scale-dependent-greta/)
    -   The code in that blog post isn't that useful, but it was an entry way to the methods from Miguet et al 2017.
-   [Miguet 2017 paper](https://doi.org/10.1111/2041-210X.12830) explaining how to estimate the decay of distance-dependent effects
-   [Moll et al. 2020](https://doi.org/10.1111/ecog.04762) on how urbanization impacts wildlife An application of the Miguet paper, including JAGS code
    -   I translated these models into Stan. I tested it with simulated data (assumed that rings of SD of NDVI predicted normally distr. response variable) and it worked well. When I tried on real data, the scale parameter wasn't informed by the data whatsoever. Maybe I need more data?
    -   I just needed more data. Everything is working well now.
-   [Monroe et al. 2022](https://doi.org/10.1002/ecs2.4320) used a scale-selecting parameter with linear interpolation to estimate the scale of effect (and uncertainty around it). Paper was on sage grouse counts and surrounding sagebrush habitat. I chatted with lead author a bit over Teams and he's readily communicative which is a plus
    -   Pretend you have a whole bunch of focal sites and you created multiple concentric buffers around each one. These buffers extend from 50 to 2000 in 50m intervals and will be used to estimate which scale of the effect is strongest. Now you have a data matrix $\mathbf{E}$ containing the mean tree mortality within each of these buffers with N rows for each focal site and M columns for each scale. In the model, you estimate a parameter $\varsigma$, bounded between 1 and M, to select which scale of $\mathbf{E}$ is strongest in predicting your response variable. Because the data is computed in discrete steps, you use linear interpolation to estimate the effect of $\mathbf{E}$. E.g. if $\varsigma$ is 3.4, its calculated as $\mathbf{E[,3]} \times .6 + \mathbf{E[,4]} \times .4$.
    -   Internally, to do this interpolation, you need to calculate the floor of $\varsigma$ and treat it as an integer in order to index $\mathbf{E}$. Stan doesn't allow parameters to be integers. I found this workaround on the Stan forum of someone else asking the exact question I needed [here](https://discourse.mc-stan.org/t/workaround-using-truncated-parameters-for-indexing/28971/2).
    -   I tested out this approach and so far it's working, even with real data. I think the approach by Miguet and Moll seems better, but this is a decent alternative if it makes the model work better.
