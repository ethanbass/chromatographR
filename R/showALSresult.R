showALSresult <- function(xals,
                          xlst,
                          tp = getTime(xals),
                          wl = getWavelength(xals),
                          mat.idx = 1:length(xlst),
                          img.col = terrain.colors(10),
                          zlim, xlab, ylab,
                          compound.col = 1:ncol(xals$S),
                          logsc = TRUE,
                          plotPureC = c("both", "spec", "conc", "none"),
                          titles, annotation = show.img,
                          PureChght = 0.33, PureCwdth = min(nplot, 5)/5 - .1,
                          show.img = TRUE) {
  def.par <- par(no.readonly = TRUE) # save default, for resetting...
  on.exit(par(def.par))
  
  plotPureC <- match.arg(plotPureC)
  
  if (missing(zlim)) {
    zlim <- range(xlst[mat.idx])
    zmin <- zlim[1]
    if (logsc) zlim <- log(zlim - zmin + 1)
  } else {
    zmin <- zlim[1]
  }
  zmin <- zlim[1]
  
  nplot <- length(mat.idx)
  
  switch(plotPureC,
         both = layout(matrix(c(1:nplot, 2*nplot + 2,
             (nplot+2):(2*nplot+1), nplot + 1),
             2, nplot+1, byrow = TRUE),
             heights = c(PureChght, 1),
             widths = rep(c(1, PureCwdth), c(nplot, 1))),
         spec = layout(matrix(c(1:nplot+1, 1), 1, nplot+1, byrow = TRUE),
             widths = rep(c(1, PureCwdth), c(nplot, 1))),
         conc = layout(matrix(c(1:nplot, nplot+1:nplot),
             2, nplot, byrow = TRUE),
             heights = c(PureChght, 1)),
         none = layout(matrix(1:nplot, nrow = 1))
         )
  
  if (!is.list(tp))
      tp <- lapply(xlst, function(x) tp)
  
  ## Titles may be given only for the plots shown, of for the
  ## whole plot sequence 
  if (missing(titles))
      titles <- names(xals$CList)
  if (length(titles) != nplot & length(titles) != length(xlst)) 
      stop("Invalid titles vector")
  if (length(titles) == nplot) {
    titles.tmp <- rep("", length(xlst))
    titles.tmp[mat.idx] <- titles
    titles <- titles.tmp
  }

  if (missing(xlab)) {
    xlab <- "Retention time (min.)"
    if (all(tp[[1]] == 1:length(tp[[1]])))
        xlab = "Retention time (idx)"
  }
  if (missing(ylab)) {
    ylab <- "Wavelength (nm)"
    if (all(wl == 1:length(wl)))
        ylab <- "Wavelength (idx)"
  }
  
  if (plotPureC == "both" | plotPureC == "conc") {
    if (all(titles == "")) {
      par(mar = c(1, 5, 2, 1), bty = "n", xaxt = "n",
          yaxs = "i", yaxt = "s")
    } else {
      par(mar = c(1, 5, 4, 1), bty = "n", xaxt = "n",
          yaxs = "i", yaxt = "s")
    }
    Cylim <- c(0, max(sapply(xals$CList[mat.idx], max)))
    for (iii in mat.idx) {
      matplot(tp[[iii]], xals$CList[[iii]], type = "l", lty = 1,
              col = compound.col, ylim = Cylim, ylab = "")
      title(main = titles[iii])
    }
  }
  
  if (plotPureC == "both" | plotPureC == "spec") {
    par(mar = c(5, 2, 1, 1), ann = FALSE, bty = "n",
        xaxs = "i", xaxt = "s", yaxt = "n")
    matplot(xals$S, wl,  type = "l", lty = 1, col = compound.col)
  }

  if (plotPureC == "both" | plotPureC == "conc" |
      all(titles == "")) { ## no titles
    par(mar = c(5, 5, 1, 1), ann = TRUE, bty = "o", xaxt = "s", yaxt = "s")
    titles <- rep("", nplot)
  } else { ## titles
    par(mar = c(5, 5, 3, 1), ann = TRUE, bty = "o", xaxt = "s", yaxt = "s")
  }
  
  for (iii in mat.idx) {
    if(logsc) {
      X <- log(xlst[[iii]] - zmin + 1)
    } else {
      X <-  xlst[[iii]]
    }
    
    if (show.img) {
      image(tp[[iii]], wl, X, col = img.col, zlim = zlim,
            xlab = xlab, ylab = ylab,
            main = titles[iii])
      box()
    } else {
      contour(tp[[iii]], wl, X, col = img.col, zlim = zlim,
              xlab = xlab, ylab = ylab,
              main = titles[iii])
      box()
    }
  }
  
  if (annotation) {
    par(mar = c(1, 1, 1, 1), ann = FALSE, bty = "n", xaxt = "n",
        yaxt = "n", xaxs = "i", yaxs = "i", xpd = NA)
    if (class(annotation) == "character") {
      plot(1:10, type = "n", axes = FALSE)
      text(5, 5, annotation, adj = c(.5, .5), cex = 1.5)
    } else {
      nBreaks <- length(img.col)
      Breaks <- seq(zlim[1], zlim[2], length = nBreaks)
      plot(NA, NA, xlim = c(0, 1), ylim = zlim * 1.1, type = "n")
      
      rect(0.9, Breaks[1:(nBreaks-1)], 1, Breaks[-1], 
           col = img.col, ##[(nBreaks-1):1], 
           if (nBreaks > 50) {border = NA})
      rect(0.9, Breaks[1], 1, Breaks[nBreaks], lwd = 1)
      par(yaxs = "i")
      
      if (logsc) {
        if (zlim[2] < 1) {
          tickmarks <- axTicks(2,
                               axp = c(zlim[1], round(zlim[2], 2), -1),
                               log = TRUE)
        } else {
          tickmarks <- axTicks(2,
                               axp = c(max(zlim[1] - 1, 1),
                                   zlim[2] - 1, 2),
                               log = TRUE)
        }
        tickmarks <- tickmarks[tickmarks > - 1]
        if (length(tickmarks) < 2)
            tickmarks <- round(zlim, 2)

        tm.pos <- log(tickmarks + 1)
      } else {
        tickmarks <- pretty(zlim)
        tm.pos <- tickmarks
      }
      
      tickmarks <- tickmarks[tm.pos < zlim[2] & tm.pos >= zlim[1]]
      tm.pos <- tm.pos[tm.pos < zlim[2] & tm.pos >= zlim[1]]
      text(.9, tm.pos, tickmarks, pos = 2)
    }
  }

  invisible()
}
