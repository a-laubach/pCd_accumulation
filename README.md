# pCd_accumulation
MATLAB implementation of the algorithm used to identify subsurface particulate Cd accumulation in the manuscript "Particulate Cadmium Accumulation in the Mesopelagic Ocean" submitted to Global Biogeochemical Cycles. 

## Description
An example implementation of an algorithm for determining the depth and magnitude of subsurface particulate element accumulation from station profiles. Accumulation is assessed as particulate concentrations in excess of those expected based on remineralization/regeneration of particles produced in the surface. Remineralization/regeneration is described using a double exponential equation fit to station profile data below the euphotic zone. Accumulation is identified via systematic removal of depths from the dataset followed by refits of the regeneration curve to remaining data. The refit with the lowest root mean square error is chosen as the apparent regeneration curve and the associated removed depths as points of accumulation. The magnitude of accumulation is defined as the offset between the measured concentration and the apparent regeneration curve at the same depth. For this implementation, cadmium-specific accumulation is confirmed by comparing to phosphorus accumulation. If phosphorus accumulation is above a defined threshold, the excess cadmium is not considered to be the result of cadmium-specific accumulation. This repository also contains particulate cadmium and phosphorus data from GEOTRACES Pacific Meridional Transect (GP15) Station 33 as an example case. 

## Getting Started

### Dependencies
This implementation requires MATLAB and the [MATLAB Curve Fitting Toolbox](https://www.mathworks.com/products/curvefitting.html).  

### Executing
By default, the implementation in `example_alg.m` loads data from `GP15_Stn33_example.mat` assuming it is in the same directory. Adjust the datapath as necessary. 

The algorithm has several calculation options that can be adjusted in the `opts` structure. By default, they are set to the values used for our calculations. Similarly, the `plt` structure has options to plot station profiles of the original fit equation and/or the final fit equation after accumulation points have been identified. By default, both options are set to 1 to create figures. Change these values to 0 to suppress figures. 

## Authors
A. Laubach: alaubach@ucsc.edu

## License
This project is licensed under the MIT License - see the LICENSE.md file for details

## Acknowledgments
GP15 Stn 33 particulate data are available as part of the [Leg 2 dataset](https://www.bco-dmo.org/dataset/919139) [[1]](#1) 

<a id="1">[1]</a>
Lam et al. (2024)
Size-fractionated major, minor, & trace particle composition and concentration from Leg 2 (Hilo, HI to Papeete, French Polynesia) of the US GEOTRACES Pacific Meridional Transect (PMT) cruise (GP15, RR1815) on R/V Roger Revelle from Oct to Nov 2018. Biological and Chemical Oceanography Data Management Office (BCO-DMO). (Version 1) Version Date 2024-01-30. doi:10.26008/1912/bco-dmo.919139.1
