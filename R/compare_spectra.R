## check peak for false 0s, etc.

compare_spectra <- function(peak, peak_table, chrom_list,
                             thresh_auto=0.95, thresh_man=NULL, r=100, plot_it=F,
                             lambda='256', zeros=F, ref='max', order_by = "distance", verbose=T,
                            plot_diff=T, ...){
  par(mfrow=c(2,1))
  if (length(ref)==1){
    ref = plot_spectrum(peak, peak_table, chrom_list, export_spectrum=T, chr=ref)
  }
  ts <- as.numeric(rownames(chrom_list[[1]]))
  lambdas <- as.numeric(colnames(chrom_list[[1]]))
  peak_tab_old <- peak_table
  if (is.null(thresh_man)){
    thresh_man <- thresh_auto
  }
  # if (length(ref>1)){
  #   ref <- sapply(ref, function(r){plot_spectrum(peak,peak_table,dat.pr,plot_chrom=F,export_spectrum=T,chr=r)[,2]})
  #   #ref1 = plot_spectrum(peak,peak_table,dat.pr,export_spectrum=T,chr=ref)[,2]
  # }
  ref.s <- rescale(ref[,1])
  #ans1 <- readline('accept peak as reference (y/n)?')
  ans1='y'
  if (ans1=='y'){
    RT <- round(peak_table['RT',peak],2)
    time <- which(elementwise.all.equal(RT,ts))+3
    if(zeros=='T'){
      chrs <- which(peak_table[,peak]==0)-3  
    } else {
      chrs <- (1:length(chrom_list))
      }
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
      spec.s <- rescale(spec)
      cor <- as.numeric(cor(ref.s, spec.s, method='pearson'))
      #cor <- cor(ref.s,spec.s,method='pearson')
      pks <- find_peaks(spec[lambda,])-1
      pks <- pks[cor[pks] > max(thresh_man,thresh_auto)]
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
        if (plot_it==T){
          matplot(c(-r:r), rescale(spec[lambda,]), type='l', 
                  xlab="", ylab='',las=2)
          matplot(c(-r:r), cor, type='l',add=T, lty=2)
          abline(v=c(-r:r)[pks], lty=3, col='lightgray')
        }
        if (cor[pk] >= thresh_auto){
          peak_table[(chr+3),peak] <- spec[lambda,pk]
          break
        }
        # else if (cor[pk]<0.5){
        #   next
        #   }
        else{
          par(mfrow=c(2,1))
          matplot(c(-r:r),rescale(spec[lambda,]),type='l',xlab='',ylab='')
          matplot(c(-r:r),cor,type='l',add=T,lty=2,col='slategray')
          abline(v=c(-r:r)[pk],lty=3,col='lightgray')
          mylabel <- substitute(corr == MYVALUE, 
                                list(MYVALUE = format(cor[pk],dig=2)))
          legend('topleft', legend = mylabel,bty='n',cex=0.6)
          legend('topright', legend = c("intensity","corr"),
                 pch = NA, lty = c(1, 3),
                 col = c(1,'slategray'), text.col = c(1,'slategray'),cex=0.4,bty='n')
          matplot(lambdas,cbind(ref.s,spec.s[,pk]),type='l',xlab='',ylab='',col=c(1,4))
          legend('topright', legend = c("reference","candidate"),
                 pch = NA, lty = c(1, 3),
                 col = c(1,4), text.col = c(1,4),cex=0.4,bty='n')
          ans2 <- readline(prompt = "accept this peak as a match (y/n)?")
          if (ans2=='y'){
            peak_table[(chr+3),peak] <- spec[lambda,pk]
            break
          } # else if (ans2=='n'){
        }
      }
    }
  }
  if (plot_diff==T){
    matplot(4:nrow(peak_table),data.frame(peak_table[-c(1:3),peak],peak_tab_old[-c(1:3),peak]), pch=20, xlab='old',ylab='new')
  }
  return(peak_table)
}
