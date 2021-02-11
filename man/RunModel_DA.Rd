\encoding{UTF-8}


\name{RunModel_DA}
\alias{RunModel_DA}
\alias{[.OutputsModelDA}
\alias{plot.OutputsModelDA}


\title{Run discharge simulations with ensemble-based data assimilation}


\usage{
RunModel_DA(InputsModel, InputsPert = NULL, Qobs = NULL,
            IndRun,
            FUN_MOD, Param,
            DaMethod = c("EnKF", "PF", "none"), NbMbr = NULL,
            StateEnKF = NULL, StatePert = NULL,
            Seed = NULL)

\method{[}{OutputsModelDA}(x, i)

\method{plot}{OutputsModelDA}(x, Qobs = NULL, main = NULL,
     ColSim = "orangered", ColObs = par("fg"), ...)
}


\arguments{
\item{InputsModel}{[list] object of class \emph{InputsModel} containing the data required to run the model  (see \code{\link[airGR]{CreateInputsModel}} for details)}

\item{InputsPert}{(optional) [list] object of class \emph{InputsModel} containing the ensembles of perturbed data required to evaluate the model probabilistic outputs, if the uncertainty in meteorological forcings is taken into account. The variables data proposed in \emph{InputsPert} will erase the variables data given in \emph{InputsModel}}

\item{Qobs}{(optional) [numeric] time series of observed flow [mm/d]}

\item{IndRun}{[numeric] index of period to be used for the model run [-]}

\item{FUN_MOD}{[function] daily hydrological model functions (GR4J, GR5J or GR6J; e.g. \code{\link[airGR]{RunModel_GR5J}}, \code{\link[airGR]{RunModel_CemaNeigeGR5J}})}

\item{Param}{[numeric] vector of model parameters (number of parameters depends on the used hydrological model, from 4 parameters in GR4J up to 10 parameters in GR6J with CemaNeige)}

\item{DaMethod}{[character] name of the data assimilation technique:
  \tabular{ll}{
    "EnKF" \tab discharge assimilation via Ensemble Kalman filter                   \cr
    "PF"   \tab discharge assimilation via Particle filter                          \cr
    "none" \tab discharge assimilation is not performed (i.e. open-loop simulation) \cr
  }}

\item{NbMbr}{(optional) [numeric] number of ensemble members (minimum of 20 recommanded for the EnKF scheme and of 30 for the PF scheme; by default=\code{NULL}: 50 if \emph{InputsPert} is not set or \emph{InputsPert$NbMbr} otherwise)}

\item{StateEnKF}{[character] vector of the names of state variables to be updated via EnKF:
  \tabular{ll}{
    "Prod" \tab level of the production store [mm]                             \cr
    "Rout" \tab level of the routing store [mm]                                \cr
    "UH1"  \tab unit hydrograph 1 levels [mm] (not defined for the GR5J model) \cr
    "UH2"  \tab unit hydrograph 2 levels [mm]
  }}

\item{StatePert}{[character] vector of the names of state variables to be perturbed via EnKF or PF:
  \tabular{ll}{
    "Prod" \tab level of the production store [mm]                             \cr
    "Rout" \tab level of the routing store [mm]                                \cr
    "UH1   \tab unit hydrograph 1 levels [mm] (not defined for the GR5J model) \cr
    "UH2"  \tab unit hydrograph 2 levels [mm]
  }}

\item{Seed}{(optional) [numeric] seed of random number generator}

\item{x}{[OutputsModelDA] containing the vector of dates (\emph{POSIXt}) and the time series of numeric values DA model outputs}

\item{i}{[integer] of the indices to subset a time series or [character] names of the elements to extract}

\item{main}{(optional) [character] an overall title for the plot (see \code{\link[graphics]{title}})}

\item{ColSim, ColObs}{(optional) [character] color to be used for simulated flow and observed flow (in any format that \code{\link[grDevices]{col2rgb}} accepts)}

\item{...}{other parameters to be passed through to plotting functions}
}


\value{
[list] runModel_DA function provides a list containing the outputs organised as follows:
  \tabular{ll}{
    \emph{$DatesR     } \tab [POSIXlt] series of dates (length of IndRun)
            \cr
            \cr
    \emph{$QsimEns    } \tab [numeric] matrix (dim(NbTime, NbMbr)) of ensemble discharge simulations
            \cr
            \cr
    \emph{$EnsStateBkg} \tab [numeric] array (dim(NbTime, NbMbr, Nstate)) of ensemble values of background model states (before the filter update)
            \cr
            \cr
    \emph{$EnsStateA  } \tab [numeric] array (dim(NbTime, NbMbr, Nstate)) of ensemble values of analysis model states (after the filter update)
            \cr
            \cr
    \emph{$NbTime     } \tab [integer] atomic vector of length of IndRun
            \cr
            \cr
    \emph{$NbMbr      } \tab [integer] atomic vector of number of ensemble members
            \cr
            \cr
    \emph{$NbState    } \tab [integer] atomic vector of number of of model states
         }

On the graphical outputs: \cr
  - solid line: medians of the input values \cr
  - polygon: minima and maxima of the input values \cr
}


\description{
Function which performs discharge ensemble simulations with the assimilation of observed discharges through the Ensemble Kalman filter (EnKF) or the Particle filter (PF) schemes. More information about the efficiency of these data assimilation schemes with GR4J can be found in Piazzi et al. (accepted).
}


\details{
Discharge observations are sequentially assimilated at each time step (if \emph{DaMethod} != "none") via the EnKF or PF schemes. Because of the sequential approach, the analysis states resulting from the assimilation procedure at time step \emph{t} are used as initial states at the following prediction time step \emph{t + 1}.
\cr\cr
When accounting for the uncertainty in model inputs (see \emph{InputsPert}), the ensemble discharge simulations are driven by perfect meteorological forecasts, which are generated according to the methodology proposed by Clark et al. (2008) (see \code{\link{CreateInputsPert}} for details).
\cr\cr
It is possible to enable/disable the update of specific state variables via EnKF by defining/not defining their names in \emph{StateEnKF}. Please note that \emph{StateEnKF} is therefore required only when using the EnKF.
\cr\cr
Both the PF and EnKF can account for uncertainty in model states by perturbing the state variables specified in \emph{StatePert} (Salamon and Feyen, 2009). It is noteworthy that, when using the EnKF, the perturbation is allowed only for the state variables included in \emph{StateEnKF}. If the state uncertainty is taken into account, the initial states at the prediction time step \emph{t + 1} are the perturbed analysis states resulting from the assimilation and perturbation procedures at time step \emph{t}.
\cr\cr
In order to ensure reproducible results, \emph{Seed} can be set to fix the randomness in the generation of perturbations.
\cr\cr
For further details and guidelines on the choice of the DA technique, see the references section.
\cr\cr
Nota: The function can be applied when using GR4J, GR5J and GR6J models (i.e. daily model time step), with or withouth the CemaNeige module (see \code{\link[airGR]{airGR}} package).
}


\examples{
library(airGRdatassim)

## loading catchment data
data(L0123001, package = "airGR")
Param <- c(X1 = 194.243, X2 = -0.088, X3 = 117.740, X4 = 1.680, X5 = 0.000)

## run period selection
IndRun <- seq(which(format(BasinObs$DatesR, format = "\%Y-\%m-\%d")=="2006-01-01"),
              which(format(BasinObs$DatesR, format = "\%Y-\%m-\%d")=="2006-01-31"))

## preparation of the InputsModel object
InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR5J, DatesR = BasinObs$DatesR,
                                 Precip = BasinObs$P, PotEvap = BasinObs$E)

## number of ensemble members
## minimum of 20 recommanded for the EnKF scheme
## minimum of 30 recommanded for the PF scheme
NbMbr <- 20L

## preparation of perturbed meteorological ensemble
InputsPert <- CreateInputsPert(FUN_MOD = RunModel_GR5J,
                               DatesR = BasinObs$DatesR,
                               Precip = BasinObs$P, PotEvap = BasinObs$E,
                               NbMbr = NbMbr)

## simulation with DA via EnKF
OutputsModelDA <- RunModel_DA(InputsModel = InputsModel,
                              InputsPert = InputsPert,
                              Qobs = BasinObs$Qmm,
                              IndRun = IndRun,
                              FUN_MOD = RunModel_GR5J, Param = Param,
                              DaMethod = "EnKF", NbMbr = NbMbr,
                              StateEnKF = c("Prod", "Rout"),
                              StatePert = c("Prod", "Rout"))

## results preview
plot(OutputsModelDA, Qobs = BasinObs$Qmm[IndRun])

## results preview on a subset
plot(OutputsModelDA[1:10], Qobs = BasinObs$Qmm[IndRun][1:10])
}


\author{
Gaia Piazzi, Olivier Delaigue
}


\references{
- Clark, M. P., Rupp, D. E., Woods, R. A., Zheng, X., Ibbitt, R. P., Slater, A. G. et al. (2008). Hydrological data assimilation with the ensemble Kalman filter: Use of streamflow observations to update states in a distributed hydrological model. Advances in Water Resources, 31(10), 1309-1324, \doi{10.1016/j.advwatres.2008.06.005}
\cr\cr
- Piazzi, G., Thirel, G., Perrin, C. and Delaigue, O. (accepted). Sequential data assimilation for streamflow forecasting: assessing the sensitivity to uncertainties and updated variables of a conceptual hydrological model. Water Resources Research, \doi{10.1029/2020WR028390}.
\cr\cr
- Salamon, P. and Feyen, L. (2009). Assessing parameter, precipitation, and predictive uncertainty in a distributed hydrological model using sequential data assimilation with the particle filter. Journal of Hydrology, 376(3-4), 428-442, \doi{10.1016/j.jhydrol.2009.07.051}
}


\seealso{
\code{\link{CreateInputsPert}}
}
