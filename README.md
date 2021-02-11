
# airGRdatassim: Suite of Tools to Perform Ensemble-Based Data Assimilation with GR Hydrological Models

## Overview

airGRdatassim is a package based on the airGR hydrological modeling package. It provides the tools to assimilate observed discharges in the GR daily hydrological model (GR4J, GR5J and GR6J, with and without the CemaNeige snow model). The package is developed at INRAE-Antony ([Catchment Hydrology research group](https://webgr.inrae.fr/en/home/) of the HYCAR Research Unit, France). More information about the efficiency of these data assimilation schemes with GR5J can be found in Piazzi et al. (accepted)


## Installation

To download the version of the airGRdatassim package that is on GitLab, you have first install the [Git software](https://git-scm.com/downloads). Since you need the latest version of airGR (not yet on CRAN), you need to install [Rtools](https://cran.r-project.org/bin/windows/Rtools/history.html) in order to do so from its sources.  Then you can install the package in the R environment, using the following command lines:

``` r
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
install.packages("remotes")
remotes::install_git(url = "https://gitlab.irstea.fr/HYCAR-Hydro/airgr")
remotes::install_git(url = "https://gitlab.irstea.fr/HYCAR-Hydro/airgrdatassim")
```


## Functions and objects

The airGRdatassim package allows users of GR Hydrological models to assimilate discharge observations with the aim of improving streamflow simulations.
The data assimilation (DA) scheme has been designed to allow the choice between two sequential ensemble-based DA techniques, namely the Ensemble Kalman filter (EnKF) and the Particle filter (PF).
The functions are coded in R and both their names and arguments are consistent with the airGR package.

With the aim of providing an user-friendly package, airGRdatassim relies on two main functions :

- `CreateInputsPerturb()` generates the probabilistic model inputs to perform the ensemble-based DA when accounting for the uncertainty in meteorological forcings;

- `RunModel_DA()` performs streamflow ensemble simulations with the assimilation of observed discharges through the EnKF or the PF scheme.

Consistently with the airGR package, both structure and class of function arguments are specifically defined to prevent the risk of mis-use and ensure the flexibility of functions. Advanced users wishing to apply the package to their own models will need to comply with these imposed structures and refer to the package source codes to get all the specification requirements.


## Hydrological model

DA schemes are designed to be coupled with GR daily hydrological model, which is implemented in the airGR package. This model can be called within the airGRdatassim package using the following airGR functions (use the command `?airGR` to get the references of the GR models): 

  - `RunModel_GR4J()`: four-parameter daily lumped hydrological model
  - `RunModel_GR5J()`: five-parameter daily lumped hydrological model
  - `RunModel_GR6J()`: six-parameter daily lumped hydrological model
  - `RunModel_CemaNeigeGR4J()`: combined use of GR4J and CemaNeige
  - `RunModel_CemaNeigeGR5J()`: combined use of GR5J and CemaNeige
  - `RunModel_CemaNeigeGR6J()`: combined use of GR6J and CemaNeige


## How to get started

Because airGRdatassim is an airGR-based package, specific airGR functions should be jointly used to ensure the proper use of the airGRdatassim tools. Indeed, before performing the DA-based streamflow simulations, the hydrological model needs to be calibrated through the airGR calibration function. Therefore, the following steps are recommended:    

  1. refer to the help for `Calibration_Michel()` in the airGR package, run the provided example and then refer to the help for `CreateCalibOptions()` to understand how a model calibration is prepared/made;
  2. refer to the help for `CreateInputsPerturb()` to understand how the probabilistic model inputs are generated, if the uncertainty in meteorological forcings is taken into account;
  3. refer to the help for `RunModel_DA()` to understand how to perform the DA-based streamflow simulations;
  4. refer to the help for `ErrorCrit_NSE()` and `CreateInputsCrit()` in the airGR package to understand how the computation of an error criterion is prepared/made.


For more information and to get started with the package, you can refer to the vignette (`vignette("get_started", package = "airGRdatassim")`).


## References

- Piazzi, G., Thirel, G., Perrin, C. and Delaigue, O. (accepted). Sequential data assimilation for streamflow forecasting: assessing the sensitivity to uncertainties and updated variables of a conceptual hydrological model. Water Resources Research, doi: [10.1029/2020WR028390](https://doi.org/10.1029/2020WR028390).
