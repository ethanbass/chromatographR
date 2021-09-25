getAllPeaks <- function (CList, lambdas, max.iter=100,...){
  if (is.numeric(lambdas)){
    lambdas <- as.character(lambdas)
  }
  peaks<-list()
  CList <- lapply(CList, function(Cmat) Cmat[,lambdas, drop=F])
  peakPositions <- lapply(CList, function(Cmat){
    apply(Cmat, 2, function(x) findpeaks(x, bounds=T))})
  result <- lapply(1:length(CList), function(smpl) {
    ptable <- lapply(1:length(peakPositions[[smpl]]), function(cmpnd){
      fitpeaks(CList[[smpl]][,cmpnd], peakPositions[[smpl]][[cmpnd]],max.iter=max.iter,...)
      })
    names(ptable) <- names(peakPositions[[smpl]])
    ptable
  })
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

## function to visually check integration accuracy
## fit is output of getAllpeaks for chrome

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
    mi <- xs[min(which(abs(diff(gaussian(xs, center=pks$rt[i], width=pks$sd[i], height = pks$height[i]))) > h*slope))]
    ma <- xs[max(which(abs(diff(gaussian(xs, center=pks$rt[i], width=pks$sd[i],height = pks$height[i]))) > h*slope))]
    if(is.na(mi)|is.na(ma)){
      next #skip bad peaks
    } else{
      xs2<-seq.int(mi,ma)
      if(fit=="gaussian"){
        polygon(new.ts[xs2], gaussian(xs2, center=mean(c(mi,ma)), width=pks$sd[i], height = pks$height[i]),
                col=scales::alpha('red',a), border=NA)
      }
      if(fit=="egh"){
        
        polygon(new.ts[xs2], egh(x=xs2, center=mean(c(mi,ma)), width=pks$sd[i], height = pks$height[i], tau=pks$tau[i]),
                col=scales::alpha('red',a), border=NA)
      }
    }
  }
}

