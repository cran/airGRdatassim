
DA_EnKF <- function(Obs, Qsim, EnsState,
                   Param, StateNames,
                   NbMbr, StateEnKF = NULL,
                   StatePert = NULL, VarThr) {

  # ------ Settings

  NbState  <- nrow(EnsState)
  IndDa    <- which(StateEnKF == 1)
  NbVarDa  <- length(IndDa)
  StateBkg <- EnsState[IndDa, , drop = FALSE]

  # member names
  MbrNames <- sprintf("Mbr_%s", seq_len(NbMbr))



  # ------ Observation error covariance matrix

  VarObs <- max(VarThr^2, (0.1*Obs)^2)

  Pert    <- rnorm(NbMbr, mean = 0, sd = sqrt(VarObs))
  ObsPert <- rep(Obs, times = NbMbr)
  ObsPert <- ObsPert + Pert
  ObsPert[ObsPert < 0] <- 0

  ObsErr <- var(Pert)



  # ------ Innovations

  Innov <- ObsPert - Qsim
  names(Innov) <- MbrNames



  # ------ Kalman Gain

  EnsMeanBkg <- rowMeans(StateBkg)
  EnsMeanQ <- mean(Qsim)

  # evaluation of anomalies
  Anom  <- as.matrix(StateBkg - EnsMeanBkg)
  AnomQ <- t(as.matrix(Qsim - EnsMeanQ))

  # evaluation of Kalman Gain
  BhtMbr  <- matrix(data = NA, nrow = NbVarDa, ncol = NbMbr,
                    dimnames = list(unlist(dimnames(StateBkg)[1], use.names = FALSE),
                                    sprintf("Mbr_%s",seq_len(NbMbr))))
  HbhtMbr <- matrix(data = NA, nrow = 1, ncol = NbMbr,
                    dimnames = list(NULL, MbrNames))

  for (iMbr in seq_len(NbMbr)) {
    BhtMbr[, iMbr] <- Anom[, iMbr] %*% t(AnomQ[iMbr])
    HbhtMbr[iMbr]  <- AnomQ[iMbr] %*% t(AnomQ[iMbr])
  }
  Bht  <- as.matrix((1/(NbMbr-1)) * rowSums(BhtMbr))
  Hbht <- as.matrix((1/(NbMbr-1)) * rowSums(HbhtMbr))

  K <- Bht %*% ((Hbht+ObsErr)^(-1))



  # ------ Analysis

  StateA <- StateBkg + K %*% Innov

  EnsStateEnkf          <- EnsState
  EnsStateEnkf[IndDa, ] <- StateA



  # ------ Constraints on analysis states

  EnsStateEnkf["Prod", EnsStateEnkf["Prod", ] < 0.05*Param[1]] <- 0.05 * Param[1]
  EnsStateEnkf["Rout", EnsStateEnkf["Rout", ] <= 0] <- 1e-3
  EnsStateEnkf["UH2" , EnsStateEnkf["UH2" , ] <  0] <- 1e-3
  if ("UH1" %in% StateNames) {
    EnsStateEnkf["UH1" , EnsStateEnkf["UH1" , ] <  0] <- 1e-3
  }

  EnsStateEnkf["Prod", EnsStateEnkf["Prod", ] > Param[1]] <- Param[1] # if Prod > X1 -> Prod = X1
  EnsStateEnkf["Rout", EnsStateEnkf["Rout", ] > Param[3]] <- Param[3] # if Rout > X3 -> Rout = X3



  # ------ States perturbation

  if (!is.null(StatePert)) {
    IndPert <- as.numeric(StateNames %in% StatePert)
    Sd0 <- apply(EnsStateEnkf, 1, sd)
    SdState <- pmin(3, pmax(1.2, Sd0))
    names(SdState) <- StateNames
    MuState  <- rep(0, times = NbState)
    names(MuState) <- StateNames
    TaoState <- matrix(data = rnorm(NbState*NbMbr, mean = MuState, sd = SdState),
                       nrow = NbState, ncol = NbMbr, byrow = TRUE,
                       dimnames = list(StateNames,
                                       MbrNames))

    TaoState[IndPert == 0, ] <- 0

    EnsStatePert <- EnsStateEnkf + TaoState

    # Positive state variables
    EnsStatePert["Prod", EnsStatePert["Prod", ] < 0.05*Param[1]] <- 0.05 * Param[1]
    EnsStatePert["Rout", EnsStatePert["Rout", ] <= 0] <- 1e-3
    EnsStatePert["UH2" , EnsStatePert["UH2" , ] <  0] <- 1e-3
    if ("UH1" %in% StateNames) {
      EnsStatePert["UH1", EnsStatePert["UH1" , ] <  0] <- 1e-3
    }

    EnsStatePert["Prod", EnsStatePert["Prod", ] > Param[1]] <- Param[1] # if Prod > X1 -> Prod = X1
    EnsStatePert["Rout", EnsStatePert["Rout", ] > Param[3]] <- Param[3] # if Rout > X3 -> Rout = X3
  }



  # ------ Outputs

  ans <- list(EnsStateEnkf = EnsStateEnkf, ObsPert = ObsPert)
  if (!is.null(StatePert)) {
    ans$EnsStatePert <- EnsStatePert
  }
  return(ans)
}
