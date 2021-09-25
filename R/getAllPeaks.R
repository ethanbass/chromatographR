
getAllPeaks <- function (CList, lambdas, max.iter=100,...){
  if (is.numeric(lambdas)){
    lambdas <- as.character(lambdas)
  }
  if (length(lambdas)>1){
  peaks<-list()
  CList2 <- lapply(CList, function(Cmat) Cmat[,lambdas])
  peakPositions <- lapply(CList2, function(Cmat){
    apply(Cmat, 2, function(x) findpeaks(x, bounds=T))})
  result <- lapply(1:length(CList2), function(smpl) {
    ptable <- lapply(1:length(peakPositions[[smpl]]), function(cmpnd){
      fitpeaks(CList2[[smpl]][,cmpnd], peakPositions[[smpl]][[cmpnd]],max.iter=max.iter,...)
      })
    names(ptable) <- names(peakPositions[[smpl]])
    ptable
  })
  } else {
    peakPositions <- lapply(CList2, function(Cmat){
      findpeaks(Cmat, bounds=T)})
    result <- lapply(1:length(CList2), function(smpl){
        fitpeaks(CList2[[smpl]], peakPositions[[smpl]],...)
      })
      #names(ptable) <- names(peakPositions[[smpl]])
      ptable
  }
  names(result) <- names(peakPositions)
  result <- lapply(result, function(smpl) lapply(smpl, function(pks) pks[apply(pks, 
                                                                               1, function(x) !any(is.na(x))), , drop = FALSE]))
  timepoints <- as.numeric(rownames(CList[[1]]))
  tdiff <- median(diff(timepoints))
  lapply(result, function(smpl) lapply(smpl, function(cmpnd) {
    x <- cmpnd
    x[, 1] <- timepoints[x[, 1]]
    x[, 2:3] <- x[, 2:3] * tdiff
    x
  }))
}

# getAllPeaks2 <- function (CList, lambdas, span = NULL, findpeaks=findpeaks_by_slope,
#                          fitpeaks=fitpeaks, ...){
#   peaks<-list()
#   CList2 <- lapply(CList, function(Cmat) Cmat[,lambdas])
#   peakPositions <- lapply(CList2, function(Cmat){
#     apply(Cmat, 2, function(x) findpeaks(x))})
#   Cmat <- CList2[5]
#   result <- lapply(1:length(CList2), function(smpl) {
#     ptable <- lapply(1:length(peakPositions[[smpl]]), function(cmpnd) fitpeaks(CList2[[smpl]], peakPositions[[smpl]][[cmpnd]],...))
#     names(ptable) <- names(peakPositions[[smpl]])
#     ptable
#   })
#   names(result) <- names(peakPositions)
#   result <- lapply(result, function(smpl) lapply(smpl, function(pks) pks[apply(pks, 
#                                                                                1, function(x) !any(is.na(x))), , drop = FALSE]))
#   timepoints <- as.numeric(rownames(CList[[1]]))
#   tdiff <- median(diff(timepoints))
#   lapply(result, function(smpl) lapply(smpl, function(cmpnd) {
#     x <- cmpnd
#     x[, 1] <- timepoints[x[, 1]]
#     x[, 2:3] <- x[, 2:3] * tdiff
#     x
#   }))
# }


## function to visually check integration accuracy
## fit is output of getallpeaks for chrome

plot_peaks <- function(chrome_list, peak_list, index=1, lambda, w=5, slope=.01,
                       points=F, a=0.5, time=c('rt','raw'),h=1){
  if (!(lambda %in% names(peak_list[[1]]))){
    stop('Error: lambda must match one of the wavelengths in your peak list')
  }
  new.ts <- as.numeric(rownames(chrome_list[[1]]))
  y<-chrome_list[[index]][,lambda]
  pks<-data.frame(peak_list[[index]][[lambda]])
  if ("tau" %in% colnames(pks)){
    fit<-"egh"
  } else{ fit<-"gaussian"}
  plot(new.ts,y,type='l',xlab='',ylab='',xaxt='n',yaxt='n')
  if (time=='rt'){
    pks$rt <- (pks$rt - new.ts[1])*100
    pks$sd <- pks$sd*100
  }
  if (points==T){
    points(new.ts[pks$rt],pks$height,pch=20,cex=0.5,col='red')
  }
  for (i in 1:nrow(pks)){
    #print(i/nrow(pks))
    xs<-seq.int((pks$rt[i]-w),(pks$rt[i]+w))
    mi <- xs[min(which(abs(diff(chromatographR:::gaussian(xs, center=pks$rt[i], width=pks$sd[i], height = pks$height[i]))) > h*slope))]
    ma <- xs[max(which(abs(diff(chromatographR:::gaussian(xs, center=pks$rt[i], width=pks$sd[i],height = pks$height[i]))) > h*slope))]
    if(is.na(mi)|is.na(ma)){
      next #skip bad peaks
    } else{
      xs2<-seq.int(mi,ma)
      if(fit=="gaussian"){
        polygon(new.ts[xs2], chromatographR:::gaussian(xs2, center=mean(c(mi,ma)), width=pks$sd[i], height = pks$height[i]),
                col=scales::alpha('red',a), border=NA)
      }
      if(fit=="egh"){
        
        polygon(new.ts[xs2], chromatographR:::egh(x=xs2, center=mean(c(mi,ma)), width=pks$sd[i], height = pks$height[i], tau=pks$tau[i]),
                col=scales::alpha('red',a), border=NA)
      }
    }
  }
}

