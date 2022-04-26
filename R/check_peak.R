#' Correct false zeros in peak table.
#' 
#' Tries to correct false zeroes for a particular peak in the
#' peak_table. This can be used as a last resort to automatically (or semi-
#' automatically) align peak data. In each chromatogram, the function compares
#' all peaks within a certain radius around the focal peak on the basis of their
#' spectral similarity to a reference spectrum.
#' 
#' @aliases check_peak
#' @importFrom graphics par matplot abline legend
#' @param peak Name of peak to be checked.
#' @param peak_table A peak_table object from \code{\link{get_peaktable}}.
#' @param chrom_list A list of chromatograms in matrix form (timepoints x
#' wavelengths).
#' @param thresh_auto Defines minimum spectral similarity threshold for
#' automatic match.
#' @param thresh_man Defines minimum spectral similarity threshold to trigger
#' manual matching process. Should be set equal or lower to `thresh_auto`.
#' @param r Defines radius around focal peak to search for matches.
#' @param plot_it Logical. If TRUE, plots spectra for comparison.
#' @param lambda Wavelength for plot.
#' @param zeros Logical. If TRUE, only check zero rows in peak table.
#' @param ref Defines chromatogram for reference spectrum. If "max" choose
#' chromatogram with maximum absorbance.
#' @param order_by If "distance", will order matches by distance from
#' theoretical retention time.
#' @param verbose Logical. If TRUE, prints verbose output to console.
#' @param plot_diff Logical.
#' @param \dots Additional arguments
#' @return A peak table similar to the input peak table, but with changes to
#' entries selected column as defined by the algorithm.
#' @author Ethan Bass
#' @seealso \code{\link{get_peaks}}
#' @noRd

check_peak <- function(peak, peak_table, chrom_list,
                        thresh_auto=0.95, thresh_man=NULL, r=100, plot_it=FALSE,
                        lambda='256', zeros=FALSE, ref='max',
                        order_by = "distance", verbose=T, plot_diff=TRUE, ...){
  if (missing(chrom_list)){
    chrom_list <- try(get(peak_table$args["chrom_list"]))
    if (inherits(chrom_list, "try-error")) stop("Chromatograms not found!")
  }
  par(mfrow=c(2,1))
  if (length(ref)==1){
    ref <- plot_spectrum(peak, peak_table, chrom_list, export_spectrum=T, chr=ref)
  }
  ts <- as.numeric(rownames(chrom_list[[1]]))
  lambdas <- as.numeric(colnames(chrom_list[[1]]))
  peak_tab_old <- peak_table
  if (is.null(thresh_man)){
    thresh_man <- thresh_auto
  }
  # if (length(ref>1)){
  #   ref <- sapply(ref, function(r){plot_spectrum(peak,peak_table,dat.pr,
  #   plot_chrom=F,export_spectrum=T,chr=r)[,2]})
  #   #ref1 = plot_spectrum(peak,peak_table,dat.pr,export_spectrum=T,chr=ref)[,2]
  # }
  ref.s <- rescale(ref[,1])
  #ans1 <- readline('accept peak as reference (y/n)?')
  ans1 <- 'y'
  if (ans1 == 'y'){
    RT <- round(peak_table['RT',peak],2)
    time <- which(elementwise.all.equal(RT,ts))+3
    if(zeros){
      chrs <- which(peak_table[,peak]==0)-3  
    } else {
      chrs <- seq_along(chrom_list)
    }
    count <- 0
    for (chr in chrs){
      count <- count+1
      ans <- 'n'
      spec <- t(chrom_list[[chr]][c((time-r):(time+r)),])
      spec.s <- rescale(spec)
      cor <- as.numeric(cor(ref.s, spec.s, method='pearson'))
      #cor <- cor(ref.s,spec.s,method='pearson')
      pks <- find_peaks(spec[lambda,])-1
      pks <- pks[cor[pks] > max(thresh_man, thresh_auto)]
      # pks <- do.call(pmax,data.frame(t(cor[,pks])))
      # do.call(which.max,data.frame(t(cor[,pks])))
      # apply(data.frame(t(cor[,pks])),1,which.max)
      if (length(pks) > 1){
        if (order_by=='height'){
          pks <- pks[order(spec[lambda,pks], decreasing=T)]}
        else if (order_by=='distance'){
          pks <- pks[order(abs(c(-r:r)[pks]))]
        }
      }
      for (pk in pks){
        if (plot_it){
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
          if (ans2 == 'y'){
            peak_table[(chr+3),peak] <- spec[lambda,pk]
            break
          } # else if (ans2=='n'){
        }
      }
    }
  }
  if (plot_diff){
    matplot(4:nrow(peak_table),data.frame(peak_table[-c(1:3),peak],peak_tab_old[-c(1:3),peak]), pch=20, xlab='old',ylab='new')
  }
  return(peak_table)
}



