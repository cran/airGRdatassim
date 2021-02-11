plot.InputsPert <- function(x, which = "all", main = NULL,
                            ColPrecip = "royalblue", ColPotEvap = "green3",
                            ask = prod(par("mfcol")) < length(which) && dev.interactive(),
                            ...) {


  ## ---------- check arguments

  ## class
  if (!inherits(x, "InputsPert")) {
    stop("'x' must be of class InputsPert")
  }

  ## which
  NamesInputsPert <- c("Precip", "PotEvap")
  which <- unique(which)
  which <- match.arg(arg = which, choices = c("all", NamesInputsPert), several.ok = TRUE)
  if (any(which %in% "all")) {
    which <- NamesInputsPert
  }
  NamesInputsPert <- intersect(names(x), which)
  if (length(NamesInputsPert) < 1L) {
    stop(sprintf("'%s' element not available in x", which))
  }

  ## ask
  if (length(NamesInputsPert) > 1L && ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }


  ## ---------- graphical variables

  TimeUnit <- c("daily", "hourly")
  TimeUnit <- match.arg(class(x), TimeUnit, several.ok = TRUE)
  TimeUnit <- switch(TimeUnit,
                     daily  = "day",
                     hourly = "hour")


  for (i in seq_along(NamesInputsPert)) {

    iName <- NamesInputsPert[i]
    IsPrecip <- iName == "Precip"
    Col  <- ifelse(test = IsPrecip, yes = ColPrecip,       no = ColPotEvap)
    YLab <- ifelse(test = IsPrecip, yes = "total precip.",  no = "pot. evap.")
    Main <- ifelse(test = IsPrecip, yes = "Precipitation", no = "Potential evapotranspiration")
    Main <- ifelse(test = is.null(main),
                   yes = sprintf("%s ensemble", Main),
                   no = main[i])


    ## ---------- plot

    dev.hold()
    RangeInputsPert <- apply(x[[iName]], MARGIN = 1, FUN = range)
    plot(x = x$DatesR, y = rowMeans(x[[iName]]),
         ylim = range(RangeInputsPert),
         type = "l", col = Col, lwd = 2,
         main = Main,
         xlab = sprintf("time [%s]", TimeUnit), ylab = sprintf("%s [mm/%s]", YLab, TimeUnit),
         panel.first = polygon(x = c(as.numeric(x$DatesR), rev(as.numeric(x$DatesR))),
                               y = c(RangeInputsPert[1L, ], rev(RangeInputsPert[2L, ])),
                               col = adjustcolor(Col, alpha.f = 0.25), border = NA))
    dev.flush()

  }


}
