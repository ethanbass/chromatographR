get_peaks <- function (chrom_list, lambdas, max.iter=100,
                         fit = c("egh", "gaussian"), sd.max=50, ...){
  fit <- match.arg(fit, c("egh", "gaussian"))
  if (is.numeric(lambdas)){
    lambdas <- as.character(lambdas)
  }
  peaks<-list()
  chrom_list <- lapply(chrom_list, function(c_mat) c_mat[,lambdas, drop=F])
  peak_positions <- lapply(chrom_list, function(c_mat){
    apply(c_mat, 2, function(x) find_peaks(x, ...))})
  result <- lapply(1:length(chrom_list), function(smpl){
    ptable <- lapply(1:length(peak_positions[[smpl]]), function(cmpnd){
      fit_peaks(chrom_list[[smpl]][,cmpnd], peak_positions[[smpl]][[cmpnd]], fit=fit, max.iter=max.iter, sd.max=sd.max)
    })
    names(ptable) <- names(peak_positions[[smpl]])
    ptable
  })
  names(result) <- names(peak_positions)
  result <- lapply(result, function(smpl) lapply(smpl, function(pks) pks[apply(pks, 
                                                                               1, function(x) !any(is.na(x))), , drop = FALSE]))
  timepoints <- as.numeric(rownames(chrom_list[[1]]))
  tdiff <- median(diff(timepoints))
  lapply(result, function(smpl) lapply(smpl, function(cmpnd) {
    x <- cmpnd
    x[, c('rt', 'start', 'end')] <- sapply(c('rt', 'start', 'end'), function(j) timepoints[x[,j]])
    x[, c('sd', 'FWHM')] <- x[, c('sd', 'FWHM')] * tdiff
    if (!is.null(x$tau)){x[, c('tau')] <- x[, c('tau')] * tdiff} 
    x
  }))
}


getAllPeaks <- function (chrom_list, lambdas, max.iter=100,
                       fit = c("egh", "gaussian"), sd.max=50, ...){
  msg<-"The function `getAllPeaks` is deprecated. Please use `get_peaks` instead"
  .Deprecated(getAllPeaks, package="chromatographR", msg,
              old = as.character(sys.call(sys.parent()))[1L])
  fit <- match.arg(fit, c("egh", "gaussian"))
  if (is.numeric(lambdas)){
    lambdas <- as.character(lambdas)
  }
  peaks<-list()
  chrom_list <- lapply(chrom_list, function(c_mat) c_mat[,lambdas, drop=F])
  peak_positions <- lapply(chrom_list, function(c_mat){
    apply(c_mat, 2, function(x) find_peaks(x, ...))})
  result <- lapply(1:length(chrom_list), function(smpl){
    ptable <- lapply(1:length(peak_positions[[smpl]]), function(cmpnd){
      fit_peaks(chrom_list[[smpl]][,cmpnd], peak_positions[[smpl]][[cmpnd]], fit=fit, max.iter=max.iter, sd.max=sd.max)
    })
    names(ptable) <- names(peak_positions[[smpl]])
    ptable
  })
  names(result) <- names(peak_positions)
  result <- lapply(result, function(smpl) lapply(smpl, function(pks) pks[apply(pks, 
                                                                               1, function(x) !any(is.na(x))), , drop = FALSE]))
  timepoints <- as.numeric(rownames(chrom_list[[1]]))
  tdiff <- median(diff(timepoints))
  lapply(result, function(smpl) lapply(smpl, function(cmpnd) {
    x <- cmpnd
    x[, c('rt', 'start', 'end')] <- sapply(c('rt', 'start', 'end'), function(j) timepoints[x[,j]])
    x[, c('sd', 'FWHM')] <- x[, c('sd', 'FWHM')] * tdiff
    if (!is.null(x$tau)){x[, c('tau')] <- x[, c('tau')] * tdiff} 
    x
  }))
}

## function to visually check integration accuracy
## fit is output of get_peaks for chrome

plot_peaks <- function(chrom_list, peak_list, index=1, lambda=NULL,
                       points=F, ticks=F, a=0.5, cex.points=0.5, ...){
  if (is.null(lambda)){
    lambda <- names(peak_list[[1]])[1]
  }
  if (!(lambda %in% names(peak_list[[1]]))){
    stop('Error: lambda must match one of the wavelengths in your peak list')
  }
  if (is.numeric(lambda)){lambda <- as.character(lambda)}
  new.ts <- as.numeric(rownames(chrom_list[[1]]))
  y <- chrom_list[[index]][,lambda]
  pks <- data.frame(peak_list[[index]][[lambda]])
  fit <- ifelse("tau" %in% colnames(pks), "egh", "gaussian")
  plot(new.ts, y, type='l', xlab='', ylab='', xaxt='n', yaxt='n', ...)
  if (points==T){
    points(pks$rt, pks$height, pch=20, cex=cex.points, col='red')
  }
  if (ticks==T){
    arrows(pks$start, y[which(new.ts %in% pks$start)]-5, pks$start, y[which(new.ts %in% pks$start)]+5, col="blue", length=0)
    arrows(pks$end, y[which(new.ts %in% pks$end)]-5, pks$end, y[which(new.ts %in% pks$end)]+5, col="blue", length=0)
  }
  for (i in 1:nrow(pks)){
    peak.loc<-seq.int((pks$start[i]),(pks$end[i]), by = .01)
      if (fit=="gaussian"){
        yvals <- gaussian(peak.loc, center=pks$rt[i], width=pks$sd[i], height = pks$height[i])
        color <- "red"
      }
      else if (fit == "egh"){
        yvals <- egh(x=peak.loc, center=pks$rt[i], width=pks$sd[i], height = pks$height[i], tau=pks$tau[i])
        color <- "purple"
      }
      sapply(1:(length(peak.loc)-1), function(i){
        polygon(peak.loc[c(i,i,(i+1),(i+1))], c(0,yvals[i:(i+1)],0),col=scales::alpha(color,a), lty=3,border=NA)
      })
  }
}
