\name{airGRdatassim}
\alias{airGRdatassim}
\docType{package}
\encoding{UTF-8}
\title{\packageTitle{airGRdatassim}}
\description{
\packageDescription{airGRdatassim} airGRdatassim currently runs on daily hydrological models (GR4J, GR5J and GR6J, with and without the CemaNeige snow model). The package is developed at INRAE-Antony (\href{https://webgr.inrae.fr/en/home/}{Catchment Hydrology research group} of the HYCAR Research Unit, France). More information about the efficiency of these data assimilation schemes with GR5J can be found in Piazzi et al. (accepted). \cr\cr


## --- Functions and objects

The airGRdatassim package allows users of GR hydrological models to assimilate discharge observations with the aim of improving streamflow simulations.
The package has been designed to allow the choice between two sequential ensemble-based data assimilation (DA) techniques, namely the Ensemble Kalman filter (EnKF) and the Particle filter (PF). \cr
The functions are coded in R and both their names and arguments are consistent with the airGR package. \cr

With the aim of providing a user-friendly package, airGRdatassim relies on two main functions: \cr
- \code{\link{CreateInputsPert}} generates the probabilistic model inputs to perform ensemble-based DA when accounting for the uncertainty in meteorological forcings; \cr
- \code{\link{RunModel_DA}} performs streamflow ensemble simulations with the assimilation of observed discharges through the EnKF or the PF schemes. \cr \cr

Consistently with the airGR package, both structure and class of function arguments are specifically defined to prevent from mis-use and to ensure the flexibility of functions. Advanced users wishing to apply the package to their own models will need to comply with these imposed structures and refer to the package source codes to get all the specification requirements. \cr\cr


## --- Models

DA schemes are designed to be coupled with GR daily hydrological models, which are implemented in the airGR package. These models can be called within the airGRdatassim package using the following airGR functions (use the command `?airGR` to get the references of the GR models): \cr
  - \code{\link[airGR]{RunModel_GR4J}}: four-parameter daily lumped hydrological model \cr
  - \code{\link[airGR]{RunModel_GR5J}}: five-parameter daily lumped hydrological model \cr
  - \code{\link[airGR]{RunModel_GR6J}}: six-parameter daily lumped hydrological model  \cr
  - \code{\link[airGR]{RunModel_CemaNeigeGR4J}}: combined use of GR4J and CemaNeige    \cr
  - \code{\link[airGR]{RunModel_CemaNeigeGR5J}}: combined use of GR5J and CemaNeige    \cr
  - \code{\link[airGR]{RunModel_CemaNeigeGR6J}}: combined use of GR6J and CemaNeige    \cr\cr


## --- How to get started

Because airGRdatassim is an airGR-based package, specific airGR functions should be jointly used to ensure the proper use of the airGRdatassim tools. Indeed, before performing the DA-based streamflow simulations, the hydrological model needs to be calibrated through the airGR calibration function. Therefore, the following steps are recommended: \cr

  1. refer to the help for \code{\link[airGR]{Calibration_Michel}} in the airGR package, run the provided example and then refer to the help for \code{\link{CreateCalibOptions}} to understand how a model calibration is prepared/made; \cr
  2. refer to the help for \code{\link{CreateInputsPert}} to understand how the probabilistic model inputs are generated, if the uncertainty in meteorological forcings is taken into account; \cr
  3. refer to the help for \code{\link{RunModel_DA}} to understand how to perform the DA-based streamflow simulations; \cr
  4. refer to the help for \code{\link[airGR]{ErrorCrit_NSE}} and \code{\link[airGR]{CreateInputsCrit}} in the airGR package to understand how the computation of an error criterion is prepared/made. \cr

For more information and to get started with the package, you can refer to the vignette: \cr
\code{vignette("get_started", package = "airGRdatassim")}. \cr\cr


## --- References

- Piazzi, G., Thirel, G., Perrin, C. and Delaigue, O. (accepted). Sequential data assimilation for streamflow forecasting: assessing the sensitivity to uncertainties and updated variables of a conceptual hydrological model. Water Resources Research, \doi{10.1029/2020WR028390}.

}



\concept{hydrology}
\concept{model}
\concept{data assimilation}
\concept{ensemble}
\concept{GR4J}
