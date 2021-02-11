## ---- warning=FALSE, include=FALSE--------------------------------------------
#knitr::write_bib("airGR", file = "airgr_ref.bib")
airGR_bib <- toBibtex(citation("airGR"))
airGR_bib <- gsub("@Article\\{", "@Article{airGR2017", airGR_bib)
airGR_bib <- gsub("@Manual\\{", "@Manual{R-airGR", airGR_bib)
airGRdatassim_bib <- toBibtex(citation("airGRdatassim")) #toBibtex(citation("airGRdatassim")[1])
airGRdatassim_bib <- gsub("@Article\\{", "@Article{piazzi_sequential_2021", airGRdatassim_bib)
airGRdatassim_bib <- gsub("@Manual\\{", "@Manual{R-airGRdatassim", airGRdatassim_bib)
options(encoding = "UTF-8")
writeLines(text = airGR_bib, con = "airgr_ref.bib")
writeLines(text = airGRdatassim_bib, con = "airgrdatassim_ref.bib")
options(encoding = "native.enc")

## ---- warning=FALSE-----------------------------------------------------------
library(airGRdatassim)

data(L0123001, package = "airGR") 
head(BasinObs)

## ---- warning=FALSE-----------------------------------------------------------
# ensemble size 
NbMbr <- 100L

# "EnKF" ; "PF" ; "none"
DaMethod <- "EnKF"

# "Prod": production store level
# "Rout": routing store level;
# "UH1": unit hydrograph 1 levels (not defined in GR5J model)
# "UH2": unit hydrograph 2 levels
StateEnKF <- c("Prod", "Rout")
StatePert <- c("Prod", "Rout") 

## ---- warning=FALSE-----------------------------------------------------------
## simulation period
IndRun <- seq(which(format(BasinObs$DatesR, format = "%Y-%m-%d") == "2006-09-01"), 
              which(format(BasinObs$DatesR, format = "%Y-%m-%d") == "2006-10-31"))

## preparation of InputsModel object
InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR5J, DatesR = BasinObs$DatesR, 
                                 Precip = BasinObs$P, PotEvap = BasinObs$E)

## generation of probabilistic model inputs
InputsPert <- CreateInputsPert(FUN_MOD = RunModel_GR5J,
                               DatesR = BasinObs$DatesR,
                               Precip = BasinObs$P, PotEvap = BasinObs$E, 
                               NbMbr = NbMbr)

## ---- warning=FALSE-----------------------------------------------------------
## discharge observations to be assimilated
Qobs <- BasinObs$Qmm[IndRun]

## ---- warning=FALSE-----------------------------------------------------------
Param <- c(X1 = 194.243, X2 = -0.088, X3 = 117.740, X4 = 1.680, X5 = 0.000)
ResEnKF <- RunModel_DA(InputsModel = InputsModel,
                       InputsPert = InputsPert,
                       Qobs = BasinObs$Qmm,
                       IndRun = IndRun,
                       FUN_MOD = RunModel_GR5J, Param = Param,
                       DaMethod = "EnKF", NbMbr = NbMbr,
                       StateEnKF = StateEnKF, StatePert = StatePert)

## ---- warning=FALSE-----------------------------------------------------------
ResPF <- RunModel_DA(InputsModel = InputsModel,
                     InputsPert = InputsPert,
                     Qobs = BasinObs$Qmm,
                     IndRun = IndRun,
                     FUN_MOD = RunModel_GR5J, Param = Param, 
                     DaMethod = "PF", NbMbr = NbMbr,
                     StatePert = "Rout")

## ---- warning=FALSE-----------------------------------------------------------
RunOptions <- CreateRunOptions(FUN_MOD = RunModel_GR5J,
                               InputsModel = InputsModel, 
                               IndPeriod_Run = IndRun)
InputsCritMulti <- CreateInputsCrit(FUN_CRIT = list(ErrorCrit_KGE, ErrorCrit_NSE), 
                                    InputsModel = InputsModel, RunOptions = RunOptions,
                                    Obs = list(Qobs, Qobs),
                                    VarObs = list("Q", "Q"), 
                                    transfo = list("", "sqrt"),
                                    Weights = NULL)

## ---- warning=FALSE, eval=TRUE, message=FALSE---------------------------------
# open-loop run (without DA)
OutputsModel_OL <- RunModel_GR5J(InputsModel = InputsModel, RunOptions = RunOptions, 
                                 Param = Param)
CritOL <- ErrorCrit(InputsCrit = InputsCritMulti, OutputsModel = OutputsModel_OL)

## ---- warning=FALSE, eval=TRUE, message=FALSE---------------------------------
# EnKF run
OutputsModel_EnKF <- OutputsModel_OL
OutputsModel_EnKF$Qsim <- rowMeans(ResEnKF$QsimEns)
CritEnKF <- ErrorCrit(InputsCrit = InputsCritMulti, OutputsModel = OutputsModel_EnKF)

## ---- warning=FALSE, eval=TRUE, message=FALSE---------------------------------
# PF run
OutputsModel_PF <- OutputsModel_OL
OutputsModel_PF$Qsim <- rowMeans(ResPF$QsimEns)
CritPF <- ErrorCrit(InputsCrit = InputsCritMulti, OutputsModel = OutputsModel_PF)

## ---- eval=TRUE, echo=FALSE---------------------------------------------------
CritVal <- rbind(
  sapply(CritOL  , "[[", "CritValue"),
  sapply(CritEnKF, "[[", "CritValue"),
  sapply(CritPF  , "[[", "CritValue")
)
CritVal <- round(CritVal, dig = 3)
colnames(CritVal) <- sapply(CritOL, "[[", "CritName")
CritVal <- cbind(data.frame(Method = c( "OL", "EnKF", "PL")),
                 CritVal)
CritVal

## ---- fig.width=7, fig.height=5, warning=FALSE, echo=FALSE--------------------
ylimMin <- min(rowMeans(ResEnKF$QsimEns), rowMeans(ResPF$QsimEns), OutputsModel_OL$Qsim, Qobs)
ylimMax <- max(rowMeans(ResEnKF$QsimEns), rowMeans(ResPF$QsimEns), OutputsModel_OL$Qsim, Qobs)

oldpar <- par(mar = c(5, 4, 2, 2))
matplot(x = BasinObs$DatesR[IndRun],
        y = cbind(OutputsModel_OL$Qsim,
                  rowMeans(ResEnKF$QsimEns),
                  rowMeans(ResPF$QsimEns),
                  Qobs),
        type = "l", col = c("red", "green", "blue", "black"),
        lty = c(1, 1, 1, 3), lwd = c(1.5, 1.5, 1.5, 1.5),
        xlab = "Time", ylab = "Discharge [mm/d]",
        xaxt = "n"
)
axis.POSIXct(side = 1, x = BasinObs$DatesR[IndRun])
legend("topleft", legend = c("OL", "EnKF", "PF", "Obs"),
       col = c("red", "green", "blue", "black"),
       lty = c(1, 1, 1, 3), lwd = rep(1.5, 4))
par(oldpar)

