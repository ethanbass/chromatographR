## check peak for false 0s, etc.

compare_spectra <- function(peak, peak_table, chrom_list, ts=new.ts, new.lambdas=new.lambdas,
                             thresh_auto=0.95, thresh_man=0.75, w=100, plot_it=F,
                             lambda='256', zeros=F, ref='max', order_by = "distance", verbose=T, ...){
  par(mfrow=c(2,1))
  if (length(ref)==1){
    ref = plot_spectrum(peak,peak_table,chrom_list,export_spectrum=T,chr=ref)[,2]
  }
  # if (length(ref>1)){
  #   ref <- sapply(ref, function(r){plot_spectrum(peak,peak_table,dat.pr,plot_chrom=F,export_spectrum=T,chr=r)[,2]})
  #   #ref1 = plot_spectrum(peak,peak_table,dat.pr,export_spectrum=T,chr=ref)[,2]
  # }
  ref.s <- scales::rescale(ref)
  ans1 <- readline('accept peak as reference (y/n)?')
  if (ans1=='y'){
    #elementwise.all.equal <- Vectorize(function(x, y) {isTRUE(all.equal(x, y))})
    RT <- round(peak_table['RT',peak],2)
    time <- which(elementwise.all.equal(RT,ts))+3
    #time <- which(ts==RT)+3
    if(zeros=='T'){
      chrs <- which(peak_table[,peak]==0)-3  
    } else {chrs <- (1:length(chrom_list))}
    #progress_bar = txtProgressBar(min=0, max=length(chrs), style = 1, char="=")
    count=0
    for (chr in chrs){
      #setTxtProgressBar(progress_bar, value = progress_bar$getVal()+1)
      count=count+1
      # if (length(chrs)%%count){
      #print(count/length(chrs))
      #svMisc::progress(chr, max.value=tail(chrs,1))
      # }
      ans='n'
      spec <- t(chrom_list[[chr]][c((time-r):(time+r)),])
      spec.s <- scales::rescale(spec)
      cor <- as.numeric(cor(ref.s,spec.s,method='pearson'))
      #cor <- cor(ref.s,spec.s,method='pearson')
      pks <- findpeaks(spec[lambda,],span=10)-1
      pks <- pks[cor[pks]>thresh_man]
      # pks <- do.call(pmax,data.frame(t(cor[,pks])))
      # do.call(which.max,data.frame(t(cor[,pks])))
      # apply(data.frame(t(cor[,pks])),1,which.max)
      if (length(pks) > 1){
        if (order_by=='height'){
          pks <- pks[order(spec[lambda,pks],decreasing=T)]}
        else if (order_by=='distance'){
          pks <- pks[order(abs(c(-r:r)[pks]))]
        }
      }
      for (pk in pks){
        # matplot(c(-r:r),scales::rescale(comp['254',]),type='l')
        # matplot(c(-r:r),c(0,scales::rescale(abs(diff(cor)))),type='l',add=T,lty=2)
        if (plot_it==T){
          matplot(c(-r:r),scales::rescale(spec[lambda,]),type='l')
          # matplot(c(-r:r),c(0,(1-scales::rescale(abs(diff(cor))))),type='l',add=T,lty=2)
          matplot(c(-r:r),cor,type='l',add=T,lty=2)
          abline(v=c(-r:r)[pks])
        }
        if (cor[pk] > thresh_auto){
          peak_table[(chr+3),peak] <- spec[lambda,pk]
          break
        }
        # else if (cor[pk]<0.5){
        #   next
        #   }
        else{
          par(mfrow=c(2,1))
          matplot(c(-r:r),scales::rescale(spec[lambda,]),type='l',xlab='',ylab='')
          # matplot(c(-r:r),c(0,(1-scales::rescale(abs(diff(cor))))),type='l',add=T,lty=2)
          matplot(c(-r:r),cor,type='l',add=T,lty=2,col='slategray')
          abline(v=c(-r:r)[pk],lty=3,col='blue')
          # abline(v=pk)
          mylabel <- substitute(corr == MYVALUE, 
                                list(MYVALUE = format(cor[pk],dig=2)))
          legend('topleft', legend = mylabel,bty='n',cex=0.6)
          legend('topright', legend = c("intensity","corr"),
                 pch = NA, lty = c(1, 3),
                 col = c(1,'slategray'), text.col = c(1,'slategray'),cex=0.4,bty='n')
          matplot(new.lambdas,cbind(ref.s,spec.s[,pk]),type='l',xlab='',ylab='',col=c(1,4))
          legend('topright', legend = c("reference","candidate"),
                 pch = NA, lty = c(1, 3),
                 col = c(1,4), text.col = c(1,4),cex=0.4,bty='n')
          ans2 <- readline(promp = "accept this peak as a match (y/n)?")
          if (ans2=='y'){
            peak_table[(chr+3),peak] <- spec[lambda,pk]
            break
          } # else if (ans2=='n'){
        }
      }
    }
  }
  return(peak_table)
}