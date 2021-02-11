
RunModel_DA <- function(InputsModel, InputsPert = NULL, Qobs = NULL,
                        IndRun,
                        FUN_MOD, Param,
                        DaMethod = c("EnKF", "PF", "none"), NbMbr = NULL,
                        StateEnKF = NULL, StatePert = NULL,
                        Seed = NULL) {

  # ------ Checks

  # FUN_MOD
  TimeUnit <- "daily"

  FUN_MODList <- c("RunModel_GR4J",
                   "RunModel_GR5J",
                   "RunModel_GR6J")

  FUN_MODSnowList <- c("RunModel_CemaNeigeGR4J",
                       "RunModel_CemaNeigeGR5J",
                       "RunModel_CemaNeigeGR6J")

  FUN_MOD <- match.fun(FUN_MOD)
  if (!any(sapply(c(FUN_MODList, FUN_MODSnowList), function(x) identical(FUN_MOD, match.fun(x))))) {
    stop(sprintf("incorrect 'FUN_MOD' for use in 'CreateInputsPerturb'. Only %s can be used",
                 paste(c(FUN_MODList, FUN_MODSnowList), collapse = ", ")))
  }

  if (identical (FUN_MOD, RunModel_GR5J) | identical (FUN_MOD, RunModel_CemaNeigeGR5J)) {
    StateNames <- c("Prod", "Rout", "UH2")
  } else {
    StateNames <- c("Prod", "Rout", "UH1", "UH2")
  }

  # DaMethod
  DaMethod <- match.arg(DaMethod)

  # Seed
  Seed0 <- as.numeric(Seed)

  # StateEnKF & StatePert
  if (DaMethod == "none" && (!is.null(StateEnKF) | !is.null(StatePert))) {
    warning("'StateEnKF' and/or 'StatePert' not taken into account when 'DaMethod' is \"none\"")
  }
  if (DaMethod == "PF" && !is.null(StateEnKF)) {
    warning("'StateEnKF' not taken into account when 'DaMethod' is \"PF\"")
  }
  if (DaMethod == "EnKF" && is.null(StateEnKF)) {
    stop("'StateEnKF' must be defined when 'DaMethod' is \"EnKF\"")
  }
  if (DaMethod != "none") {
    if (!is.null(StateEnKF)) {
      StateEnKF <- match.arg(StateEnKF, choices = StateNames, several.ok = TRUE)
    }
    if (!is.null(StatePert)) {
      StatePert <- match.arg(StatePert, choices = StateNames, several.ok = TRUE)
    }
  }
  if (DaMethod == "EnKF" && any(!StatePert %in% StateEnKF)) {
    stop(sprintf("Perturbation is allowed only for the state variables updated via EnKF (%s). Please check the consistency between 'StatePert' and 'StateEnKF'",
                 sQuote(paste(StateEnKF, collapse = ", "))))
  }
  if (DaMethod != "none" && any(!StatePert %in% StateNames)) {
    warning(StatePert[!StatePert %in% StateNames])
  }

  # InputsModel
  if (!inherits(InputsModel, "InputsModel")) {
    stop("'InputsModel' must of class 'InputsModel'")
  }

  # InputsPert
  if (is.null(InputsPert)) {
    IsMeteo <- FALSE
    if (is.null(NbMbr)) {
      NbMbr <- 50
      message("'InputsPert' and 'NbMbr' not defined: number of ensemble members automatically set to 50")
    }
  } else {
    if (!inherits(InputsPert, "InputsPert")) {
      stop("'InputsPert' must of class 'InputsPert' or NULL")
    } else {
      if (length(InputsPert$DatesR) != length(InputsModel$DatesR)) {
        stop("'InputsPert' elements must have the same length as the 'InputsModel'elements")
      }
      IsMeteo <- TRUE
      ClassInputsModel <- class(InputsModel)
      ClassInputsPert  <- class(InputsPert)[-grep("InputsPert", class(InputsPert))]
      ClassDiffInputsModel <- setdiff(ClassInputsModel, ClassInputsPert)
      ClassDiffInputsPert  <- setdiff(ClassInputsPert, ClassInputsModel)
      if (length(ClassDiffInputsModel) != 0 | length(ClassDiffInputsPert) != 0) {
        msgClassInputs <- "'InputsModel' and 'InputsPert' classes are not consistent:"
        if (length(ClassDiffInputsModel) != 0) {
          msgClassInputs <- sprintf("%s\n\tInputsModel: %s", msgClassInputs, paste(dQuote(ClassDiffInputsModel), collapse = "\t"))
        }
        if (length(ClassDiffInputsPert) != 0) {
          msgClassInputs <- sprintf("%s\n\tInputsPert:  %s", msgClassInputs, paste(dQuote(ClassDiffInputsPert), collapse = "\t"))
        }
        stop(msgClassInputs)
      }
    }

    # NbMbr
    if (is.null(NbMbr)) {
      message(sprintf("'NbMbr' not defined: number of ensemble members automatically set to 'InputsPert$NbMbr' (%i)", InputsPert$NbMbr))
      NbMbr <- InputsPert$NbMbr
    }
    if (!(is.atomic(NbMbr) && is.numeric(NbMbr) && length(NbMbr) == 1 && NbMbr >= 2)) {
      stop("'NbMbr' should be a single vector of a numeric value >= 2")
    } else {
      NbMbr <- as.integer(NbMbr)
      if (IsMeteo) {
        NbMbrMeteo <- ncol(InputsPert[[2]])
        if (NbMbr > NbMbrMeteo) {
          stop(sprintf("cannot take a number of ensemble members (%i) larger than the number available for the perturbed meteorological variables (%i)",
                       NbMbr, NbMbrMeteo))
        }
        if (NbMbr < NbMbrMeteo) {
          warning(sprintf("only %i ensemble members are taken, whereas the number available for the perturbed meteorological variables is equal to %i",
                          NbMbr, NbMbrMeteo))
        }
      }
    }

    # Qobs
    if (is.null(Qobs) || all(is.na(Qobs)) || all(Qobs < 0,  na.rm = TRUE)) {
      DaMethod <- "none"
      warning("'DaMethod' is automatically set to 'none'. All Qobs may be 'NULL', 'NA' or negative")
    } else {
      if (length(Qobs) != length(InputsModel$DatesR)) {
        stop("'Qobs' must have the same length as the 'InputsModel' elements")
      }
      if (any(Qobs < 0,  na.rm = TRUE)) {
        warning("negative value(s) of Qobs are automatically set to 'NA'")
      }
    }

  }


  # ------ Settings

  # data assimilation method used (not open-loop simulation)
  IsDa <- DaMethod != "none"

  NbTime <- length(IndRun)

  NbState <- length(StateNames)

  Qobs[Qobs < 0] <- NaN
  VarThr <- quantile(Qobs, probs = 0.1, na.rm = TRUE)

  # member names
  MbrNames <- sprintf("Mbr_%s", seq_len(NbMbr))

  # time names
  TimeNames <- sprintf("Time_%s", seq_len(NbTime))

  # InputsModel
  InputsModel <- InputsModel[IndRun]

  # InputsPert
  InputsPert <- InputsPert[IndRun]

  # Qobs
  Qobs <- Qobs[IndRun]



  # ------ Ensemble initializations

  ObsPert <- matrix(data = NA,
                    nrow = NbMbr, ncol = NbTime,
                    dimnames = list(MbrNames,
                                    TimeNames))

  QsimEns <- ObsPert

  IniStatesEns   <- list()
  IniStatesEnsNbTime <- list()

  EnsStateBkg <- array(data = rep(NaN, times = NbState*NbMbr*NbTime),
                       dim = c(NbState, NbMbr, NbTime),
                       dimnames = list(StateNames,
                                       MbrNames,
                                       TimeNames))

  EnsStateA <- EnsStateBkg

  ItAssim <- 0

  # fake RunOptions
  RunOptionsIni <- airGR::CreateRunOptions(FUN_MOD = FUN_MOD,
                                           InputsModel = InputsModel,
                                           IndPeriod_Run = 1L,
                                           warning = FALSE, verbose = FALSE)
  RunOptionsIter <- airGR::CreateRunOptions(FUN_MOD = FUN_MOD,
                                            InputsModel = InputsModel,
                                            IndPeriod_Run = 1L,
                                            IndPeriod_WarmUp = 0L,
                                            IniStates = NULL,
                                            warning = FALSE, verbose = FALSE)



  # ------ Run

  if (IsMeteo) {
    if (is.null(InputsPert$Precip)) {
      InputsPert$Precip <- replicate(n = ncol(InputsPert$PotEvap),
                                     expr = InputsModel$Precip)
      dimnames(InputsPert$Precip) <- dimnames(InputsPert$PotEvap)
    }
    if (is.null(InputsPert$PotEvap)) {
      InputsPert$PotEvap <- replicate(n = ncol(InputsPert$Precip),
                                      expr = InputsModel$PotEvap)
      dimnames(InputsPert$PotEvap) <- dimnames(InputsPert$Precip)
    }
  }

  for (iTime in seq_along(IndRun)) {
    if (!is.null(Seed)) {
      Seed <- Seed0 + iTime
      set.seed(Seed)
      on.exit(set.seed(seed = NULL))
    }
    for (iMbr in seq_len(NbMbr)) {

      if (iTime == 1) { # default (one year by default) warmup

        RunOptionsIni$IndPeriod_Run <- iTime
        OutputsModel <- FUN_MOD(InputsModel = InputsModel,
                                RunOptions = RunOptionsIni, Param = Param)

      } else { # IF iTime > 1

        IniStates <- IniStatesEns[[iMbr]]
        IniStates$Store$Rest <- rep(NA, times = 3)
        IniStates <- unlist(IniStates)
        IniStates[is.na(IniStates)] <- 0
        RunOptionsIter$IniStates <- IniStates
        RunOptionsIter$IniResLevels <- NULL

        # definition of run options
        if (IsMeteo) {
          InputsPertMbr <- InputsPert
          InputsPertMbr$Precip <- InputsPert$Precip[, iMbr]
          InputsPertMbr$PotEvap <- InputsPert$PotEvap[, iMbr]
          RunOptionsIter$IndPeriod_Run <- as.integer(iTime)
          InputsModel <- InputsPertMbr
        } else {
          RunOptionsIter$IndPeriod_Run <- as.integer(iTime)
        } # END IF(IsMeteo)
        OutputsModel <- FUN_MOD(InputsModel = InputsModel,
                                RunOptions = RunOptionsIter, Param = Param)
      } # END IF(t == 1)

      IniStatesEns[[iMbr]] <- OutputsModel$StateEnd
      names(IniStatesEns)[iMbr] <- sprintf("Mbr_%s", iMbr)

      EnsStateBkg["Prod", iMbr, iTime] <- OutputsModel$Prod
      EnsStateBkg["Rout", iMbr, iTime] <- OutputsModel$Rout
      EnsStateBkg["UH2" , iMbr, iTime] <- OutputsModel$StateEnd$UH$UH2[1]
      if ("UH1" %in% StateNames) {
        EnsStateBkg["UH1", iMbr, iTime] <- OutputsModel$StateEnd$UH$UH1[1]
      }

      QsimEns[iMbr, iTime]  <- OutputsModel$Qsim

    } # END FOR particles



    # ------ Assimilation [if an observation is available]

    if (IsDa & is.finite(Qobs[iTime])) {

      ItAssim <- ItAssim + 1

      if (DaMethod == "EnKF") {
        ans <- DA_EnKF(Obs = Qobs[iTime], Qsim = QsimEns[, iTime], EnsState = EnsStateBkg[, , iTime],
                       Param = Param, StateNames = StateNames,
                       StatePert = StatePert,
                       NbMbr = NbMbr,
                       StateEnKF = StateEnKF, VarThr = VarThr)

        for (iMbr in seq_len(NbMbr)) { # olivier, it is possible to write the following 3 loops without loops?
          IniStatesEns[[iMbr]]$Store$Prod <- ans$EnsStateEnkf["Prod", iMbr]
          IniStatesEns[[iMbr]]$Store$Rout <- ans$EnsStateEnkf["Rout", iMbr]
          IniStatesEns[[iMbr]]$UH$UH2[1]  <- ans$EnsStateEnkf["UH2" , iMbr]
          if ("UH1" %in% StateNames) {
            IniStatesEns[[iMbr]]$UH$UH1[1] <- ans$EnsStateEnkf["UH1", iMbr]
          }
        }

        if (iTime < NbTime) {
          IniStatesEnsNbTime[[iTime+1]] <- IniStatesEns
          names(IniStatesEnsNbTime)[iTime+1] <- sprintf("Time_%s",iTime+1)
        }

        if (!is.null(StatePert)) {
          for (iMbr in seq_len(NbMbr)) {
            IniStatesEns[[iMbr]]$Store$Prod <- ans$EnsStatePert["Prod", iMbr]
            IniStatesEns[[iMbr]]$Store$Rout <- ans$EnsStatePert["Rout", iMbr]
            IniStatesEns[[iMbr]]$UH$UH2[1]  <- ans$EnsStatePert["UH2" , iMbr]
            if ("UH1" %in% StateNames) {
              IniStatesEns[[iMbr]]$UH$UH1[1] <- ans$EnsStatePert["UH1", iMbr]
            }
          }
        }

        EnsStateA[, , iTime] <- ans$EnsStateEnkf

        if (iTime < NbTime) {
          if (!is.null(StatePert)) {
            EnsStateBkg[, , iTime+1] <- ans$EnsStatePert
          } else {
            EnsStateBkg[, , iTime+1] <- ans$EnsStateEnkf
          }
        }
        ObsPert[, iTime] <- ans$ObsPert

      } else if (DaMethod == "PF") {
        ans <- DA_PF(Obs = Qobs[iTime], Qsim = QsimEns[, iTime], States = IniStatesEns,
                     Param = Param, StateNames = StateNames,
                     NbMbr = NbMbr,
                     StatePert = StatePert, VarThr = VarThr)

        if (!is.null(StatePert)) {
          IniStatesEns <- ans$EnsStatePert
        } else {
          IniStatesEns <- ans$EnsStatePf
        }

        EnsStateA["Prod", , iTime] <- sapply(seq_along(ans$EnsStatePf), function(x) ans$EnsStatePf[[x]]$Store$Prod)
        EnsStateA["Rout", , iTime] <- sapply(seq_along(ans$EnsStatePf), function(x) ans$EnsStatePf[[x]]$Store$Rout)
        EnsStateA["UH2" , , iTime] <- sapply(seq_along(ans$EnsStatePf), function(x) ans$EnsStatePf[[x]]$UH$UH2[1])
        if ("UH1" %in% StateNames) {
          EnsStateA["UH1", , iTime] <- sapply(seq_along(ans$EnsStatePf), function(x) ans$EnsStatePf[[x]]$UH$UH1[1])
        }

        if (iTime < NbTime) { # olivier?
          IniStatesEnsNbTime[[iTime+1]] <- ans$EnsStatePf
          names(IniStatesEnsNbTime)[iTime+1] <- sprintf("Time_%s", iTime+1)
        }
      }

    } else { # IF no assimilation

      if (iTime < NbTime) {
        IniStatesEnsNbTime[[iTime+1]] <- IniStatesEns
        names(IniStatesEnsNbTime)[iTime+1] <- sprintf("Time_%s", iTime+1)

        EnsStateBkg[, , iTime+1] <- EnsStateBkg[, , iTime]
      }

      EnsStateA [, , iTime] <- EnsStateBkg[, , iTime]
      if (DaMethod == "EnKF") {
        ObsPert[, iTime] <- rep(Qobs[iTime], times = NbMbr)
      }

    } # END IF assimilation

  } # END FOR time



  # ------ Outputs and class

  res <- list(DatesR = InputsModel$DatesR,
              QsimEns = t(QsimEns),
              EnsStateBkg = aperm(EnsStateBkg),
              EnsStateA = aperm(EnsStateA),
              NbTime = NbTime,
              NbMbr = NbMbr,
              NbState = NbState)
  class(res) <- c("OutputsModelDA", "OutputsModel", DaMethod, TimeUnit)
  return(res)

}
