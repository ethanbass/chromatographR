get_peaks <- function (chrom_list, lambdas, max.iter=100,
                         fit = c("egh", "gaussian", "emg"), ...){
  fit <- match.arg(fit, c("egh", "gaussian", "emg"))
  if (is.numeric(lambdas)){
    lambdas <- as.character(lambdas)
  }
  peaks<-list()
  chrom_list <- lapply(chrom_list, function(c_mat) c_mat[,lambdas, drop=F])
  peak_positions <- lapply(chrom_list, function(c_mat){
    apply(c_mat, 2, function(x) find_peaks(x, bounds=T))})
  result <- lapply(1:length(chrom_list), function(smpl) {
    ptable <- lapply(1:length(peak_positions[[smpl]]), function(cmpnd){
      fit_peaks(chrom_list[[smpl]][,cmpnd], peak_positions[[smpl]][[cmpnd]], fit=fit, max.iter=max.iter, ...)
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
                         fit = c("egh", "gaussian", "emg"), ...){
  msg<-"The function `getAllPeaks` is deprecated. Please use `get_peaks` instead"
  .Deprecated(getAllPeaks, package="chromatographR", msg,
              old = as.character(sys.call(sys.parent()))[1L])
  fit <- match.arg(fit, c("egh", "gaussian", "emg"))
  if (is.numeric(lambdas)){
    lambdas <- as.character(lambdas)
  }
  peaks<-list()
  chrom_list <- lapply(chrom_list, function(c_mat) c_mat[,lambdas, drop=F])
  peak_positions <- lapply(chrom_list, function(c_mat){
    apply(c_mat, 2, function(x) find_peaks(x, bounds=T))})
  result <- lapply(1:length(chrom_list), function(smpl) {
    ptable <- lapply(1:length(peak_positions[[smpl]]), function(cmpnd){
      fit_peaks(chrom_list[[smpl]][,cmpnd], peak_positions[[smpl]][[cmpnd]], fit=fit, max.iter=max.iter, ...)
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
## fit is output of getAllpeaks for chrome

plot_peaks <- function(chrom_list, peak_list, index=1, lambda=NULL, w=0.5, slope=.01, h=1,
                       points=F, a=0.5){
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
  if ("tau" %in% colnames(pks)){
    fit <- "egh"
  } else{ fit <- "gaussian"}
  plot(new.ts, y, type='l', xlab='', ylab='', xaxt='n', yaxt='n')
  if (points==T){
    points(pks$rt, pks$height, pch=20, cex=0.5, col='red')
  }
  for (i in 1:nrow(pks)){
    xs<-seq.int((pks$rt[i]-w),(pks$rt[i]+w), by = .01)
    m <- gaussian(xs, center=pks$rt[i], width=pks$sd[i], height = pks$height[i])
    mi <- xs[min(which(abs(diff(m)) > h*slope))]
    ma <- xs[max(which(abs(diff(m)) > h*slope))]
    if(is.na(mi)|is.na(ma)){
      next #skip bad peaks
    } else{
      xs2<-seq.int(mi,ma, by = .01)
      if(fit=="gaussian"){
        polygon(xs2, gaussian(xs2, center=pks$rt[i], width=pks$sd[i], height = pks$height[i]),
                col=scales::alpha('red',a), border=NA)
      }
      if(fit=="egh"){
        polygon(xs2, egh(x=xs2, center=pks$rt[i], width=pks$sd[i], height = pks$height[i], tau=pks$tau[i]),
                col=scales::alpha('purple',a), border=NA)
      }
    }
  }
}

