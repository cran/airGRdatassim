
## =================================================================================
## function to extract parts of InputsPert or OutputsModelDA objects
## =================================================================================

'[.InputsPert' <-  function(x, i) {
  res <- NextMethod()
  if (is.numeric(i)) {
    res$NbMbr <- x$NbMbr
  }
  return(res)
}


'[.OutputsModelDA' <- function(x, i) {
  res <- NextMethod()
  if (is.numeric(i)) {
    res <- airGR:::.ExtractOutputsModel(x, i)
    res$NbMbr   <- x$NbMbr
    res$NbTime  <- x$NbTime
    res$NbState <- x$NbState
  }
  return(res)
}
